// Assembles the single-page app from src/web into app/index.html (desktop, offline fonts)
// and dist/stackwise-web.html (browser version, Google Fonts). Run: npm run build-web
import { readFileSync, writeFileSync, mkdirSync } from 'node:fs';
const r = p => readFileSync(new URL('../src/web/' + p, import.meta.url), 'utf8').replace(/\r\n/g, '\n');
const tpl = r('app.template.html');
const head = tpl.slice(0, tpl.indexOf('</style>'));
let page = head + r('extra.css') + '</style>\n\n' + r('body.html');
for (const [k, f] of [['PARAMS', 'params.json'], ['ENGINE', 'engine.js'], ['STRUCT', 'struct.js'], ['APP', 'app.js']]) page = page.replace(`/*${k}*/`, () => r(f).trim());
const base = '<style>[hidden]{display:none!important}*{box-sizing:border-box}body{margin:0}</style>';
const gf = /<link rel="preconnect"[^>]*>\n<link rel="stylesheet" href="https:\/\/fonts\.googleapis\.com[^>]*>\n/;
if (!gf.test(page)) throw new Error('font link not found');
const wrap = (body, extraHead) => `<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n${extraHead}\n${base}\n</head>\n<body>\n${body}\n</body>\n</html>\n`;
const desktop = wrap(page.replace(gf, ''), '<link rel="stylesheet" href="fonts.css">');
writeFileSync(new URL('../app/index.html', import.meta.url), desktop);
mkdirSync(new URL('../dist/', import.meta.url), { recursive: true });
writeFileSync(new URL('../dist/stackwise-web.html', import.meta.url), wrap(page, ''));
console.log('app/index.html', desktop.length, 'bytes');
