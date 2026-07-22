// Regenerate the homepage tile screenshots in img/. Needs playwright-core
// (npm i playwright-core) and a local server: python3 -m http.server 8642.
// Then: node scripts/capture_screenshots.js && sips -Z 1000 -s format jpeg img/*.png ...
const { chromium } = require('playwright-core');
const os = require('os');
const path = require('path');

const EXEC = path.join(os.homedir(),
  'Library/Caches/ms-playwright/chromium_headless_shell-1208/chrome-headless-shell-mac-arm64/chrome-headless-shell');

const BASE = 'http://localhost:8642';
const OUT = '/Users/fields/discgolf/img';
const SHOTS = [
  { name: 'designer',    url: BASE + '/disc_designer.html', wait: 1500 },
  { name: 'wind_tunnel', url: BASE + '/wind_tunnel.html',   wait: 9000 },
  { name: 'database',    url: BASE + '/database.html',      wait: 1500 },
  { name: 'flight_sim',  url: BASE + '/flight_sim.html',    wait: 1500 },
  { name: 'docs',        url: BASE + '/docs.html',          wait: 2000 },
];

(async () => {
  const browser = await chromium.launch({ executablePath: EXEC });
  const page = await browser.newPage({ viewport: { width: 1280, height: 800 }, deviceScaleFactor: 2 });
  for (const s of SHOTS) {
    await page.goto('about:blank');
    await page.goto(s.url, { waitUntil: 'networkidle' });
    await page.waitForTimeout(s.wait);
    await page.screenshot({ path: `${OUT}/${s.name}.png` });
    console.log('shot', s.name);
  }
  await browser.close();
})();
