# Stackwise

Oligo thermodynamics bench for Mac and Windows: double-stranded sequence maps, nearest-neighbor Tm, hairpins, self- and cross-dimers, and primer checks. Works fully offline.

## What it does

- **Sequences:** paste or open a strand (plain, FASTA, GenBank). The complement is drawn underneath in 10-bp blocks with position numbers. Drag across either strand to select a fragment, see its Tm, GC, ΔG, hairpin and self-dimer live, analyze it without saving, trim or extend its ends one base at a time, then save it.
- **Oligos:** Tm (SantaLucia 1998 unified nearest-neighbor parameters, SantaLucia or Owczarzy salt correction), ΔH/ΔS/ΔG, melting curve, MW, ε260, hairpins, self- and cross-dimers drawn as duplexes, 3′ stability profile and primer checks.
- **Conditions:** Na⁺, Mg²⁺, dNTP, oligo concentration and ΔG temperature apply everywhere.

Tm matches Primer3 to within 0.01 °C. Hairpin and dimer ΔG use Primer3's nearest-neighbor parameter tables with a dynamic-programming search and agree with Primer3 to about 0.3–0.6 kcal/mol on average.

## Layout

| Path | What it is |
| --- | --- |
| `src/web/` | Source of the app: page shell, styles, thermodynamics engine (`engine.js`), structure search (`struct.js`), UI (`app.js`), parameter tables (`params.json`) |
| `scripts/build-web.mjs` | Assembles `src/web` into `app/index.html` (desktop, bundled fonts) and `dist/stackwise-web.html` (browser) |
| `app/` | What the desktop app loads: `index.html`, `fonts.css`, `fonts/` |
| `src-tauri/` | Tauri shell: window, bundle and icon settings |
| `app-icon.png` | Source for all icons (`npx tauri icon app-icon.png -o src-tauri/icons`) |

Your library is stored inside the app (it stays between launches).

## Running it locally

Needs Node 18+ and Rust (https://rustup.rs). On a Mac, also the Xcode command line tools (`xcode-select --install`).

```sh
npm install
npm run build-web   # after editing anything in src/web
npm run dev         # opens the app
npm run build       # makes the installer for the computer you run it on
```

## Builds

GitHub builds the installers on every push to `main`. Look under **Actions → Build Stackwise → latest run → Artifacts**:

- `Stackwise-mac` contains a `.dmg` that runs on both Apple Silicon and Intel Macs.
- `Stackwise-windows` contains a `.exe` installer and an `.msi`.

Pushing a tag such as `v0.1.0` also creates a draft Release with the installers attached.

### First launch (the app isn't signed)

- **Mac:** open the `.dmg` and drag Stackwise into Applications. The first time you open it, macOS blocks it. Go to **System Settings → Privacy & Security** and click **Open Anyway**. If macOS instead says the app "is damaged", run `xattr -cr /Applications/Stackwise.app` in Terminal once.
- **Windows:** if SmartScreen shows "Windows protected your PC", click **More info → Run anyway**.
