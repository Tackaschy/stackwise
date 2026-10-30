// Stackwise is a self-contained web app (app/index.html). The native shell hosts it and
// adds one thing a web page cannot do: talk to NCBI, whose servers do not accept
// cross-origin requests from pages. Only the endpoints listed below can be reached.

fn endpoint(name: &str) -> Option<&'static str> {
    match name {
        "blast" => Some("https://blast.ncbi.nlm.nih.gov/Blast.cgi"),
        "blastalign" => Some("https://blast.ncbi.nlm.nih.gov/BlastAlign.cgi"),
        "primerblast" => Some("https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi"),
        _ => None,
    }
}

/// Sends one request to NCBI and returns the response body.
/// `method` is "get" (query string), "post" (urlencoded form) or "multipart" (form-data).
#[tauri::command]
async fn ncbi_request(endpoint_name: String, method: String, params: Vec<(String, String)>) -> Result<String, String> {
    let url = endpoint(&endpoint_name).ok_or("Unknown NCBI service")?;
    let client = reqwest::Client::builder()
        .user_agent(concat!("Stackwise/", env!("CARGO_PKG_VERSION")))
        .timeout(std::time::Duration::from_secs(180))
        .build()
        .map_err(|e| e.to_string())?;
    let req = match method.as_str() {
        "get" => client.get(url).query(&params),
        "post" => client.post(url).form(&params),
        "multipart" => {
            let mut form = reqwest::multipart::Form::new();
            for (k, v) in params {
                form = form.text(k, v);
            }
            client.post(url).multipart(form)
        }
        _ => return Err("Unknown request method".into()),
    };
    let resp = req.send().await.map_err(|e| format!("Could not reach NCBI: {e}"))?;
    let status = resp.status();
    let body = resp.text().await.map_err(|e| e.to_string())?;
    if !status.is_success() {
        return Err(format!("NCBI answered {status}"));
    }
    Ok(body)
}

/// Kept for the BLAST pane of 0.2.0: CMD=Put is posted, everything else is a GET.
#[tauri::command]
async fn ncbi_blast(params: Vec<(String, String)>) -> Result<String, String> {
    let is_put = params.iter().any(|(k, v)| k == "CMD" && v == "Put");
    ncbi_request("blast".into(), if is_put { "post" } else { "get" }.into(), params).await
}

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    tauri::Builder::default()
        .plugin(tauri_plugin_opener::init())
        .invoke_handler(tauri::generate_handler![ncbi_request, ncbi_blast])
        .run(tauri::generate_context!())
        .expect("error while running Stackwise");
}
