// Stackwise is a self-contained web app (app/index.html). The native shell hosts it and
// adds one thing a web page cannot do: talk to NCBI BLAST, which does not allow
// cross-origin requests from pages.

const BLAST_URL: &str = "https://blast.ncbi.nlm.nih.gov/Blast.cgi";

/// Sends one request to the NCBI BLAST Common URL API and returns the response body.
/// `CMD=Put` is sent as a form POST (queries can be long); everything else as a GET.
#[tauri::command]
async fn ncbi_blast(params: Vec<(String, String)>) -> Result<String, String> {
    let client = reqwest::Client::builder()
        .user_agent(concat!("Stackwise/", env!("CARGO_PKG_VERSION")))
        .timeout(std::time::Duration::from_secs(120))
        .build()
        .map_err(|e| e.to_string())?;
    let is_put = params.iter().any(|(k, v)| k == "CMD" && v == "Put");
    let req = if is_put {
        client.post(BLAST_URL).form(&params)
    } else {
        client.get(BLAST_URL).query(&params)
    };
    let resp = req.send().await.map_err(|e| format!("Could not reach NCBI: {e}"))?;
    let status = resp.status();
    let body = resp.text().await.map_err(|e| e.to_string())?;
    if !status.is_success() {
        return Err(format!("NCBI answered {status}"));
    }
    Ok(body)
}

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    tauri::Builder::default()
        .plugin(tauri_plugin_opener::init())
        .invoke_handler(tauri::generate_handler![ncbi_blast])
        .run(tauri::generate_context!())
        .expect("error while running Stackwise");
}
