'use strict';
/* Meteor trajectory viewer: finds the newest KML/KMZ in the repo, fits the trajectory, and builds a combined KML
   (trajectory + visibility circles + sight lines/observer pin) for the phone's location.  Port of meteor_trajectory.ipynb. */

const REPO = 'JRRowe/Impactor';
const CFG = {
  fitTolKm: 0.5,        // max distance from the fitted line for a point to count as on the trajectory
  minAltKm: 10,         // points below this are ignored (observers, dropper feet)
  refraction: 0.13,     // standard refraction coefficient (effective Earth radius = R / (1 - k))
  nCircle: 360,
  lineWidth: 6, markerScale: 1.5, pointIconScale: 1.0,
};
// WGS84
const A = 6378137.0, F = 1 / 298.257223563, E2 = F * (2 - F), B = A * (1 - F);
const rad = d => d * Math.PI / 180, deg = r => r * 180 / Math.PI;

// ---------- vector helpers ----------
const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
const mul = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const norm = a => Math.hypot(a[0], a[1], a[2]);

// ---------- geodesy ----------
function geodeticToEcef(lat, lon, alt) {
  const p = rad(lat), l = rad(lon), n = A / Math.sqrt(1 - E2 * Math.sin(p) ** 2);
  return [(n + alt) * Math.cos(p) * Math.cos(l), (n + alt) * Math.cos(p) * Math.sin(l), (n * (1 - E2) + alt) * Math.sin(p)];
}
function ecefToGeodetic(v) {
  const [x, y, z] = v, lon = Math.atan2(y, x), p = Math.hypot(x, y);
  let lat = Math.atan2(z, p * (1 - E2)), n, alt;
  for (let i = 0; i < 10; i++) {
    n = A / Math.sqrt(1 - E2 * Math.sin(lat) ** 2);
    alt = p * Math.cos(lat) + z * Math.sin(lat) - A * A / n;
    lat = Math.atan2(z, p * (1 - E2 * n / (n + alt)));
  }
  n = A / Math.sqrt(1 - E2 * Math.sin(lat) ** 2);
  alt = p * Math.cos(lat) + z * Math.sin(lat) - A * A / n;
  return { lat: deg(lat), lon: deg(lon), alt };
}
function lookAngles(lat0, lon0, alt0, target) {
  const d = sub(target, geodeticToEcef(lat0, lon0, alt0));
  const la = rad(lat0), lo = rad(lon0);
  const east = [-Math.sin(lo), Math.cos(lo), 0];
  const north = [-Math.sin(la) * Math.cos(lo), -Math.sin(la) * Math.sin(lo), Math.cos(la)];
  const up = [Math.cos(la) * Math.cos(lo), Math.cos(la) * Math.sin(lo), Math.sin(la)];
  const e = dot(d, east), n = dot(d, north), u = dot(d, up);
  return { az: (deg(Math.atan2(e, n)) + 360) % 360, el: deg(Math.atan2(u, Math.hypot(e, n))), range: norm(d) };
}
function vincentyDirect(lat, lon, azDeg, dist) {
  const al1 = rad(azDeg), tanU1 = (1 - F) * Math.tan(rad(lat));
  const cosU1 = 1 / Math.sqrt(1 + tanU1 * tanU1), sinU1 = tanU1 * cosU1;
  const sig1 = Math.atan2(tanU1, Math.cos(al1));
  const sinAl = cosU1 * Math.sin(al1), cos2Al = 1 - sinAl * sinAl;
  const u2 = cos2Al * (A * A - B * B) / (B * B);
  const AA = 1 + u2 / 16384 * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)));
  const BB = u2 / 1024 * (256 + u2 * (-128 + u2 * (74 - 47 * u2)));
  let sig = dist / (B * AA);
  for (let i = 0; i < 100; i++) {
    const c2sm = Math.cos(2 * sig1 + sig), s = Math.sin(sig), c = Math.cos(sig);
    const dsig = BB * s * (c2sm + BB / 4 * (c * (-1 + 2 * c2sm * c2sm) - BB / 6 * c2sm * (-3 + 4 * s * s) * (-3 + 4 * c2sm * c2sm)));
    const nw = dist / (B * AA) + dsig, done = Math.abs(nw - sig) < 1e-12;
    sig = nw;
    if (done) break;
  }
  const c2sm = Math.cos(2 * sig1 + sig), s = Math.sin(sig), c = Math.cos(sig);
  const tmp = sinU1 * s - cosU1 * c * Math.cos(al1);
  const lat2 = Math.atan2(sinU1 * c + cosU1 * s * Math.cos(al1), (1 - F) * Math.sqrt(sinAl * sinAl + tmp * tmp));
  const lam = Math.atan2(s * Math.sin(al1), cosU1 * c - sinU1 * s * Math.cos(al1));
  const C = F / 16 * cos2Al * (4 + F * (4 - 3 * cos2Al));
  const L = lam - (1 - C) * F * sinAl * (sig + C * s * (c2sm + C * c * (-1 + 2 * c2sm * c2sm)));
  return [deg(lat2), ((lon + deg(L) + 540) % 360 + 360) % 360 - 180];
}
function vincentyInverse(lat1, lon1, lat2, lon2) {   // -> [distance m, initial azimuth deg]
  if (lat1 === lat2 && lon1 === lon2) return [0, 0];
  const L = rad(lon2 - lon1);
  const u1 = Math.atan((1 - F) * Math.tan(rad(lat1))), u2_ = Math.atan((1 - F) * Math.tan(rad(lat2)));
  const su1 = Math.sin(u1), cu1 = Math.cos(u1), su2 = Math.sin(u2_), cu2 = Math.cos(u2_);
  let lam = L, ss, cs, sig, sa, c2a, c2sm, converged = false;
  for (let i = 0; i < 200; i++) {
    const sl = Math.sin(lam), cl = Math.cos(lam);
    ss = Math.sqrt((cu2 * sl) ** 2 + (cu1 * su2 - su1 * cu2 * cl) ** 2);
    if (ss === 0) return [0, 0];
    cs = su1 * su2 + cu1 * cu2 * cl;
    sig = Math.atan2(ss, cs);
    sa = cu1 * cu2 * sl / ss; c2a = 1 - sa * sa;
    c2sm = c2a !== 0 ? cs - 2 * su1 * su2 / c2a : 0;
    const C = F / 16 * c2a * (4 + F * (4 - 3 * c2a)), prev = lam;
    lam = L + (1 - C) * F * sa * (sig + C * ss * (c2sm + C * cs * (-1 + 2 * c2sm * c2sm)));
    if (Math.abs(lam - prev) < 1e-12) { converged = true; break; }
  }
  if (!converged) {                                  // near-antipodal: spherical fallback (far outside any circle)
    const p1 = rad(lat1), p2 = rad(lat2), dl = rad(lon2 - lon1);
    const h = Math.sin((p2 - p1) / 2) ** 2 + Math.cos(p1) * Math.cos(p2) * Math.sin(dl / 2) ** 2;
    const brg = Math.atan2(Math.sin(dl) * Math.cos(p2), Math.cos(p1) * Math.sin(p2) - Math.sin(p1) * Math.cos(p2) * Math.cos(dl));
    return [2 * 6371008.8 * Math.asin(Math.min(1, Math.sqrt(h))), (deg(brg) + 360) % 360];
  }
  const u2 = c2a * (A * A - B * B) / (B * B);
  const AA = 1 + u2 / 16384 * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)));
  const BB = u2 / 1024 * (256 + u2 * (-128 + u2 * (74 - 47 * u2)));
  const dsig = BB * ss * (c2sm + BB / 4 * (cs * (-1 + 2 * c2sm * c2sm) - BB / 6 * c2sm * (-3 + 4 * ss * ss) * (-3 + 4 * c2sm * c2sm)));
  const az1 = Math.atan2(cu2 * Math.sin(lam), cu1 * su2 - su1 * cu2 * Math.cos(lam));
  return [B * AA * (sig - dsig), (deg(az1) + 360) % 360];
}
function horizonDistance(lat, alt, azDeg) {          // metres along the ellipsoid, with standard refraction
  const phi = rad(lat), al = rad(azDeg), w = Math.sqrt(1 - E2 * Math.sin(phi) ** 2);
  const nR = A / w, mR = A * (1 - E2) / w ** 3;
  const rAz = 1 / (Math.cos(al) ** 2 / mR + Math.sin(al) ** 2 / nR);
  const re = rAz / (1 - CFG.refraction);
  return re * Math.acos(re / (re + alt));
}

// ---------- reading KML / KMZ ----------
async function unzipKml(buf) {
  const dv = new DataView(buf), u8 = new Uint8Array(buf);
  let e = buf.byteLength - 22;
  while (e >= 0 && dv.getUint32(e, true) !== 0x06054b50) e--;
  if (e < 0) throw new Error('KMZ is not a valid zip file');
  const n = dv.getUint16(e + 10, true);
  let p = dv.getUint32(e + 16, true);
  const entries = [];
  for (let i = 0; i < n && dv.getUint32(p, true) === 0x02014b50; i++) {
    const method = dv.getUint16(p + 10, true), csize = dv.getUint32(p + 20, true);
    const nlen = dv.getUint16(p + 28, true), xlen = dv.getUint16(p + 30, true), clen = dv.getUint16(p + 32, true);
    const off = dv.getUint32(p + 42, true);
    entries.push({ name: new TextDecoder().decode(u8.subarray(p + 46, p + 46 + nlen)), method, csize, off });
    p += 46 + nlen + xlen + clen;
  }
  const ent = entries.find(x => /(^|\/)doc\.kml$/i.test(x.name)) || entries.find(x => /\.kml$/i.test(x.name));
  if (!ent) throw new Error('No .kml inside the KMZ');
  const start = ent.off + 30 + dv.getUint16(ent.off + 26, true) + dv.getUint16(ent.off + 28, true);
  const data = u8.subarray(start, start + ent.csize);
  if (ent.method === 0) return new TextDecoder().decode(data);
  const stream = new Blob([data]).stream().pipeThrough(new DecompressionStream('deflate-raw'));
  return new TextDecoder().decode(await new Response(stream).arrayBuffer());
}

function parseKml(text) {
  const doc = new DOMParser().parseFromString(text, 'application/xml');
  if (doc.getElementsByTagName('parsererror').length) throw new Error('Could not parse the KML');
  const q = (el, name) => Array.from(el.getElementsByTagNameNS('*', name));
  const verts = [];
  for (const pm of q(doc, 'Placemark')) {
    let name = '', when = null;
    for (const ch of pm.children) {
      if (ch.localName === 'name') name = ch.textContent.trim();
      else if (ch.localName === 'TimeStamp') { const w = q(ch, 'when')[0]; if (w) when = w.textContent.trim(); }
    }
    for (const kind of ['Point', 'LineString', 'LinearRing']) {
      for (const g of q(pm, kind)) {
        const c = q(g, 'coordinates')[0];
        if (!c) continue;
        const am = q(g, 'altitudeMode')[0], mode = am ? am.textContent.trim() : 'clampToGround';
        for (const tok of c.textContent.trim().split(/\s+/)) {
          const v = tok.split(',').map(Number);
          if (v.length < 2 || isNaN(v[0]) || isNaN(v[1])) continue;
          verts.push({ lon: v[0], lat: v[1], alt: v[2] || 0, name, when, kind, mode });
        }
      }
    }
    for (const g of q(pm, 'Track')) {
      const whens = q(g, 'when').map(w => w.textContent.trim());
      q(g, 'coord').forEach((c, k) => {
        const v = c.textContent.trim().split(/\s+/).map(Number);
        if (v.length >= 2) verts.push({ lon: v[0], lat: v[1], alt: v[2] || 0, name, when: whens[k] || when, kind: 'Track', mode: 'absolute' });
      });
    }
  }
  return verts;
}

// ---------- trajectory fit ----------
function mulberry32(a) { return () => { a |= 0; a = a + 0x6D2B79F5 | 0; let t = Math.imul(a ^ a >>> 15, 1 | a); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; }; }
function distToLine(p, p0, u) { return norm(cross(sub(p, p0), u)); }
function principalDir(pts, c) {                       // dominant eigenvector of the covariance (power iteration)
  const C = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (const p of pts) { const d = sub(p, c); for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) C[i][j] += d[i] * d[j]; }
  let v = sub(pts[pts.length - 1], pts[0]);
  if (norm(v) === 0) v = [1, 0, 0];
  v = mul(v, 1 / norm(v));
  for (let it = 0; it < 100; it++) {
    const w = [dot(C[0], v), dot(C[1], v), dot(C[2], v)], nw = norm(w);
    if (nw === 0) break;
    v = mul(w, 1 / nw);
  }
  return v;
}
function meanOf(pts) { let s = [0, 0, 0]; for (const p of pts) s = add(s, p); return mul(s, 1 / pts.length); }

function fitTrajectory(P, tolM) {
  const n = P.length;
  if (n < 2) throw new Error('Need at least two points above the altitude cut-off to define a trajectory.');
  const pairs = [];
  if (n * (n - 1) / 2 <= 20000) { for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) pairs.push([i, j]); }
  else { const r = mulberry32(1); for (let k = 0; k < 20000; k++) pairs.push([Math.floor(r() * n), Math.floor(r() * n)]); }
  let best = null, bestCnt = -1, bestRms = Infinity;
  for (const [i, j] of pairs) {
    if (i === j) continue;
    const d = sub(P[j], P[i]), L = norm(d);
    if (L < 2000) continue;
    const u = mul(d, 1 / L);
    let cnt = 0, ss = 0; const mask = new Array(n);
    for (let k = 0; k < n; k++) { const dist = distToLine(P[k], P[i], u); mask[k] = dist <= tolM; if (mask[k]) { cnt++; ss += dist * dist; } }
    const rms = Math.sqrt(ss / Math.max(cnt, 1));
    if (cnt > bestCnt || (cnt === bestCnt && rms < bestRms)) { best = mask; bestCnt = cnt; bestRms = rms; }
  }
  let mask = best || P.map(() => true);
  const sel = m => P.filter((_, k) => m[k]);
  for (let it = 0; it < 10; it++) {
    const pts = sel(mask), c = meanOf(pts), u = principalDir(pts, c);
    const nm = P.map(p => distToLine(p, c, u) <= tolM);
    if (nm.filter(Boolean).length < 2 || nm.every((v, k) => v === mask[k])) break;
    mask = nm;
  }
  const pts = sel(mask), c = meanOf(pts), u = principalDir(pts, c);
  const ts = pts.map(p => dot(sub(p, c), u));
  let ends = [add(c, mul(u, Math.min(...ts))), add(c, mul(u, Math.max(...ts)))];
  if (ecefToGeodetic(ends[0]).alt < ecefToGeodetic(ends[1]).alt) ends.reverse();
  const resid = pts.map(p => distToLine(p, c, u));
  return { mask, start: ends[0], end: ends[1], mid: mul(add(ends[0], ends[1]), 0.5),
           rmsM: Math.sqrt(resid.reduce((a, b) => a + b * b, 0) / resid.length), maxM: Math.max(...resid) };
}

// ---------- KML writers ----------
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
function kStyle(id, o = {}) {
  const { line, width = 3, fill, icon, scale = 1, labels = true } = o;
  let s = `<Style id="${id}">`;
  if (!labels) s += '<LabelStyle><scale>0</scale></LabelStyle>';
  if (icon) s += `<IconStyle><color>${icon}</color><scale>${scale}</scale><Icon><href>http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png</href></Icon></IconStyle>`;
  if (line) s += `<LineStyle><color>${line}</color><width>${width}</width></LineStyle>`;
  if (fill !== undefined) s += `<PolyStyle><color>${fill}</color></PolyStyle>`;
  return s + '</Style>\n';
}
function pmPoint(name, lon, lat, alt, o = {}) {
  const { desc = '', style = '', mode = 'absolute', when = null, extrude = false } = o;
  let s = '<Placemark>' + (name ? `<name>${esc(name)}</name>` : '');
  if (desc) s += `<description>${esc(desc)}</description>`;
  if (when) s += `<TimeStamp><when>${esc(when)}</when></TimeStamp>`;
  if (style) s += `<styleUrl>#${style}</styleUrl>`;
  return s + `<Point>${extrude ? '<extrude>1</extrude>' : ''}<altitudeMode>${mode}</altitudeMode><coordinates>${lon.toFixed(8)},${lat.toFixed(8)},${alt.toFixed(2)}</coordinates></Point></Placemark>\n`;
}
function pmLine(name, pts, o = {}) {                  // pts: [[lat, lon, alt], ...]
  const { desc = '', style = '', mode = 'absolute', extrude = false } = o;
  const coords = pts.map(([la, lo, al]) => `${lo.toFixed(8)},${la.toFixed(8)},${al.toFixed(2)}`).join(' ');
  let s = `<Placemark><name>${esc(name)}</name>`;
  if (desc) s += `<description>${esc(desc)}</description>`;
  if (style) s += `<styleUrl>#${style}</styleUrl>`;
  return s + `<LineString>${extrude ? '<extrude>1</extrude>' : ''}<tessellate>0</tessellate><altitudeMode>${mode}</altitudeMode><coordinates>${coords}</coordinates></LineString></Placemark>\n`;
}
function pmPolygon(name, ring, o = {}) {              // ring: [[lat, lon], ...]
  const { desc = '', style = '' } = o;
  const coords = ring.map(([la, lo]) => `${lo.toFixed(6)},${la.toFixed(6)},0`).join(' ');
  let s = `<Placemark><name>${esc(name)}</name>`;
  if (desc) s += `<description>${esc(desc)}</description>`;
  if (style) s += `<styleUrl>#${style}</styleUrl>`;
  return s + `<Polygon><tessellate>1</tessellate><altitudeMode>clampToGround</altitudeMode><outerBoundaryIs><LinearRing><coordinates>${coords}</coordinates></LinearRing></outerBoundaryIs></Polygon></Placemark>\n`;
}
const kmlDoc = (name, body) => `<?xml version="1.0" encoding="UTF-8"?>\n<kml xmlns="http://www.opengis.net/kml/2.2">\n<Document>\n<name>${esc(name)}</name>\n${body}</Document>\n</kml>\n`;
const folder = (name, content) => `<Folder><name>${esc(name)}</name>\n${content}</Folder>\n`;

// ---------- event (fixed part: trajectory + circles) ----------
function circleRing(lat, lon, alt) {
  const ring = []; let sum = 0;
  for (let k = 0; k < CFG.nCircle; k++) {
    const az = 360 * k / CFG.nCircle, r = horizonDistance(lat, alt, az);
    sum += r; ring.push(vincentyDirect(lat, lon, az, r));
  }
  ring.push(ring[0]);
  return { ring, meanKm: sum / CFG.nCircle / 1000 };
}

function computeEvent(verts) {
  const cand = verts.filter(v => v.alt >= CFG.minAltKm * 1000);
  const P = cand.map(v => geodeticToEcef(v.lat, v.lon, v.alt));
  const fit = fitTrajectory(P, CFG.fitTolKm * 1000);
  const used = cand.filter((_, k) => fit.mask[k]);
  const S = ecefToGeodetic(fit.start), E = ecefToGeodetic(fit.end), M = ecefToGeodetic(fit.mid);
  const lengthKm = norm(sub(fit.end, fit.start)) / 1000;
  const top = circleRing(S.lat, S.lon, S.alt), bot = circleRing(E.lat, E.lon, E.alt);

  const nseg = Math.max(2, Math.ceil(lengthKm) + 1), linePts = [];
  for (let k = 0; k < nseg; k++) { const g = ecefToGeodetic(add(fit.start, mul(sub(fit.end, fit.start), k / (nseg - 1)))); linePts.push([g.lat, g.lon, g.alt]); }
  let tb = kStyle('traj', { line: 'ff0000ff', width: CFG.lineWidth }) + kStyle('pt', { icon: 'ff00ffff' }) +
           kStyle('key', { icon: 'ff0000ff', scale: CFG.markerScale });
  tb += pmLine('Computed trajectory', linePts, { desc: `Length ${lengthKm.toFixed(1)} km; start ${(S.alt / 1000).toFixed(1)} km, end ${(E.alt / 1000).toFixed(1)} km altitude`, style: 'traj', extrude: true });
  tb += pmLine('Trajectory ground track', linePts, { style: 'traj', mode: 'clampToGround' });
  tb += pmPoint('Start (highest point)', S.lon, S.lat, S.alt, { style: 'key', extrude: true });
  tb += pmPoint('Midpoint', M.lon, M.lat, M.alt, { style: 'key', extrude: true });
  tb += pmPoint('End (lowest point)', E.lon, E.lat, E.alt, { style: 'key', extrude: true });
  tb += '<Folder><name>Input points used in the fit</name>\n';
  for (const v of used) tb += pmPoint(v.name, v.lon, v.lat, v.alt, { style: 'pt', when: v.when, extrude: true });
  tb += '</Folder>\n';

  let cb = kStyle('big', { line: 'ff00ffff', width: 3, fill: '2200ffff' }) + kStyle('small', { line: 'ffffff00', width: 3, fill: '22ffff00' }) +
           kStyle('ctr', { icon: 'ffffffff' });
  cb += pmPolygon('Visibility of HIGHEST point', top.ring, { desc: `Mean radius ${top.meanKm.toFixed(1)} km (point at ${(S.alt / 1000).toFixed(1)} km altitude)`, style: 'big' });
  cb += pmPolygon('Visibility of LOWEST point', bot.ring, { desc: `Mean radius ${bot.meanKm.toFixed(1)} km (point at ${(E.alt / 1000).toFixed(1)} km altitude)`, style: 'small' });
  cb += pmPoint('Ground point below highest point', S.lon, S.lat, 0, { mode: 'clampToGround', style: 'ctr' });
  cb += pmPoint('Ground point below lowest point', E.lon, E.lat, 0, { mode: 'clampToGround', style: 'ctr' });

  return { fit, S, E, M, lengthKm, nCand: cand.length, nUsed: used.length, top, bot, trajBody: tb, circBody: cb };
}

// ---------- observer (phone-specific part) ----------
function computeObserver(ev, lat, lon, alt) {
  const inside = (c, calt) => { const [d, az] = vincentyInverse(c.lat, c.lon, lat, lon); return d <= horizonDistance(c.lat, calt, az); };
  const inTop = inside(ev.S, ev.S.alt), inBot = inside(ev.E, ev.E.alt);
  const topIsBig = ev.S.alt >= ev.E.alt;
  const inBig = topIsBig ? inTop : inBot, inSmall = topIsBig ? inBot : inTop;
  const out = { inBig, inSmall, rows: [], pointsKml: null };
  let combined = folder('Computed trajectory', ev.trajBody) + folder('Visibility circles', ev.circBody);

  if (!inBig) {
    combined += folder('Observer', kkObs() + pmPoint('Observer (trajectory not visible from here)', lon, lat, alt,
      { mode: 'clampToGround', style: 'obs', desc: `Latitude ${lat}, longitude ${lon}: outside the larger visibility circle` }));
  } else {
    const targets = [['Bottom (end)', ev.fit.end], ['Midpoint', ev.fit.mid], ['Top (start)', ev.fit.start]];
    let sb = kStyle('sl', { line: 'ff00ff00', width: CFG.lineWidth }) + kStyle('obs', { icon: 'ff00ff00', scale: CFG.markerScale }) +
             kStyle('tgt', { icon: 'ff00ff00', scale: CFG.markerScale });
    sb += pmPoint('Observer', lon, lat, alt, { style: 'obs' });
    for (const [label, ecef] of targets) {
      const la = lookAngles(lat, lon, alt, ecef), g = ecefToGeodetic(ecef);
      const desc = `Azimuth ${la.az.toFixed(2)} deg, altitude ${la.el.toFixed(2)} deg, length ${(la.range / 1000).toFixed(2)} km`;
      out.rows.push({ label, az: la.az, el: la.el, km: la.range / 1000 });
      sb += pmLine(`To ${label}`, [[lat, lon, alt], [g.lat, g.lon, g.alt]], { desc, style: 'sl' });
      sb += pmPoint(label, g.lon, g.lat, g.alt, { style: 'tgt', extrude: true, desc });
    }
    combined += folder('Sight lines from observer', sb);
    out.pointsKml = sightlinePoints(targets, lat, lon, alt);
  }
  out.combinedKml = kmlDoc('Combined output - meteor trajectory, visibility and sight lines', combined);
  return out;
}
function kkObs() { return kStyle('obs', { icon: 'ff00ff00' }); }

function sightlinePoints(targets, lat, lon, alt) {
  const o = geodeticToEcef(lat, lon, alt), dists = [];
  for (let k = 1; k <= 90; k++) dists.push(k / 3);
  for (let k = 31; k <= 100; k++) dists.push(k);
  for (let k = 105; k <= 300; k += 5) dists.push(k);
  const colours = { 'Bottom (end)': 'ff0000ff', 'Midpoint': 'ff00ffff', 'Top (start)': 'ffffff00' };
  let body = kStyle('obs', { icon: 'ff00ff00' }) + pmPoint('Observer', lon, lat, alt, { style: 'obs' });
  for (const [label, ecef] of targets) {
    const sid = 'p' + label.replace(/\W/g, '');
    body += kStyle(sid, { icon: colours[label], scale: CFG.pointIconScale, labels: false }) + kStyle(sid + 'L', { line: colours[label], width: 2 });
    const d = sub(ecef, o), u = mul(d, 1 / norm(d));
    let pts = '';
    for (const x of dists) { const g = ecefToGeodetic(add(o, mul(u, x))); pts += pmPoint('', g.lon, g.lat, g.alt, { style: sid }); }
    body += folder(`Points along line to ${label}`, pts);
    const end = ecefToGeodetic(add(o, mul(u, dists[dists.length - 1])));
    body += pmLine(`Guide line to ${label}`, [[lat, lon, alt], [end.lat, end.lon, end.alt]], { style: sid + 'L' });
  }
  return kmlDoc('Points along sight lines', body);
}

// ---------- finding the newest file in the repo ----------
async function findNewestFile() {
  const params = new URLSearchParams(location.search);
  if (params.get('url')) return { name: params.get('url').split('/').pop().split('?')[0], url: params.get('url') };
  const res = await fetch(`https://api.github.com/repos/${REPO}/contents/`);
  if (!res.ok) throw new Error(`GitHub listing failed (${res.status}). Try again later.`);
  const files = (await res.json()).filter(f => f.type === 'file' && /\.(kml|kmz)$/i.test(f.name));
  if (!files.length) throw new Error('No KML or KMZ files found in the repository.');
  // newest by the timestamp in the file name (YYYYMMDDThhmmss...); otherwise by last commit date
  for (const f of files) { const m = f.name.match(/^(\d{8}T\d{6})/); f.key = m ? m[1] : null; }
  for (const f of files.filter(f => !f.key)) {
    try {
      const r = await fetch(`https://api.github.com/repos/${REPO}/commits?path=${encodeURIComponent(f.path)}&per_page=1`);
      const j = await r.json();
      f.key = (j[0].commit.committer.date || '').replace(/[-:]/g, '');
    } catch (e) { f.key = ''; }
  }
  files.sort((a, b) => (a.key < b.key ? 1 : a.key > b.key ? -1 : 0));
  return { name: files[0].name, url: files[0].download_url };
}

async function loadEvent() {
  const f = await findNewestFile();
  const res = await fetch(f.url);
  if (!res.ok) throw new Error(`Could not download ${f.name} (${res.status}).`);
  const text = /\.kmz$/i.test(f.name) ? await unzipKml(await res.arrayBuffer()) : await res.text();
  const ev = computeEvent(parseKml(text));
  ev.fileName = f.name;
  return ev;
}

// ---------- page ----------
const $ = id => document.getElementById(id);
let EVENT = null, lastUrls = [];

function parseAngle(s, limit) {
  const m = /^\s*([+-]?\d+(?:\.\d*)?)\s*°?\s*([NSEWnsew])?\s*$/.exec(s || '');
  if (!m) return null;
  let v = parseFloat(m[1]);
  if (m[2] && 'SWsw'.includes(m[2])) v = -Math.abs(v);
  return Math.abs(v) <= limit ? v : null;
}
async function lookupElevation(lat, lon) {
  const tries = [
    [`https://api.open-meteo.com/v1/elevation?latitude=${lat.toFixed(5)}&longitude=${lon.toFixed(5)}`, j => j.elevation[0]],
    [`https://api.opentopodata.org/v1/srtm30m?locations=${lat.toFixed(5)},${lon.toFixed(5)}`, j => j.results[0].elevation],
  ];
  for (const [url, pick] of tries) {
    try { const r = await fetch(url); if (r.ok) { const v = pick(await r.json()); if (v !== null && v !== undefined) return Number(v); } } catch (e) { /* next */ }
  }
  return null;
}
function setLink(a, kml, name) {
  const blob = new Blob([kml], { type: 'application/vnd.google-earth.kml+xml' }), url = URL.createObjectURL(blob);
  lastUrls.push(url); a.href = url; a.download = name;
}

async function run(lat, lon, altIn, gpsAlt) {
  $('status').textContent = 'Calculating...'; $('status').className = 'small';
  try {
    let alt = altIn, altSrc = 'typed';
    if (alt === null) { alt = await lookupElevation(lat, lon); altSrc = 'Open-Meteo elevation database'; }
    if (alt === null && gpsAlt !== null && gpsAlt !== undefined) { alt = gpsAlt; altSrc = 'phone GPS'; }
    if (alt === null) { alt = 0; altSrc = 'assumed 0 m (lookup failed)'; }
    const r = computeObserver(EVENT, lat, lon, alt);
    const stem = EVENT.fileName.replace(/\.(kml|kmz)$/i, '');
    let h = `<b>Location</b> ${lat.toFixed(5)}, ${lon.toFixed(5)}; height ${alt.toFixed(0)} m (${altSrc}).<br>`;
    if (!r.inBig) {
      h += '<b>trajectory is not visible from location</b><br>Location: latitude ' + lat + ', longitude ' + lon +
           '.<br>The KML shows the trajectory, the visibility circles and a pin at your location.';
    } else {
      h += 'Inside the larger visibility circle.' + (r.inSmall ? '' : ' (Outside the smaller circle: the lowest part is below your horizon.)');
      h += '<table><tr><th>Line to</th><th>Azimuth&deg;</th><th>Altitude&deg;</th><th>Length km</th></tr>' +
           r.rows.map(x => `<tr><td>${x.label}</td><td>${x.az.toFixed(2)}</td><td>${x.el.toFixed(2)}</td><td>${x.km.toFixed(2)}</td></tr>`).join('') + '</table>';
    }
    $('summary').innerHTML = h;
    lastUrls.forEach(URL.revokeObjectURL); lastUrls = [];
    const fname = `${stem}_combined_output.kml`;
    setLink($('dl'), r.combinedKml, fname);
    $('dlpts').classList.toggle('hidden', !r.pointsKml);
    if (r.pointsKml) setLink($('dlpts'), r.pointsKml, `${stem}_sightline_points.kml`);
    $('share').onclick = async () => {
      const file = new File([r.combinedKml], fname, { type: 'application/vnd.google-earth.kml+xml' });
      try {
        if (navigator.canShare && navigator.canShare({ files: [file] })) { await navigator.share({ files: [file], title: fname }); return; }
      } catch (e) { if (e.name === 'AbortError') return; }
      $('dl').click();
    };
    $('result').classList.remove('hidden');
    $('status').textContent = '';
  } catch (e) { $('status').textContent = 'Error: ' + e.message; $('status').className = 'small err'; }
}

$('go').onclick = () => {
  if (!navigator.geolocation) { $('status').textContent = 'This browser cannot give its location. Enter it by hand below.'; return; }
  $('status').textContent = 'Getting your location...';
  navigator.geolocation.getCurrentPosition(
    p => run(p.coords.latitude, p.coords.longitude, null, p.coords.altitude),
    e => { $('status').textContent = 'Could not get your location (' + e.message + '). Enter it by hand below.'; $('status').className = 'small err'; },
    { enableHighAccuracy: true, timeout: 20000, maximumAge: 0 });
};
$('manual').onclick = () => {
  const lat = parseAngle($('lat').value, 90), lon = parseAngle($('lon').value, 180);
  if (lat === null || lon === null) { $('status').textContent = 'Please enter a valid latitude (-90..90) and longitude (-180..180).'; $('status').className = 'small err'; return; }
  const a = $('alt').value.trim();
  run(lat, lon, a === '' ? null : parseFloat(a), null);
};

loadEvent().then(ev => {
  EVENT = ev;
  $('event').innerHTML = `<b>${esc(ev.fileName)}</b><br>Trajectory: ${ev.nUsed} of ${ev.nCand} points used (fit RMS ${ev.fit.rmsM.toFixed(0)} m).<br>` +
    `Start ${(ev.S.alt / 1000).toFixed(1)} km &rarr; end ${(ev.E.alt / 1000).toFixed(1)} km altitude, length ${ev.lengthKm.toFixed(1)} km.<br>` +
    `Visible out to about ${ev.top.meanKm.toFixed(0)} km (top) / ${ev.bot.meanKm.toFixed(0)} km (bottom) from the ground track.`;
  $('go').disabled = false; $('manual').disabled = false;
}).catch(e => { $('event').innerHTML = `<span class="err">Could not load the trajectory: ${esc(e.message)}</span>`; });
