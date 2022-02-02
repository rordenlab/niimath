const fs = require('fs');

let fileConents = fs.readFileSync('package.json');
let jdata = JSON.parse(fileConents);
newVersion = `0.0.${Date.now()}`
jdata['version'] = newVersion
let outdata = JSON.stringify(jdata, null, 2);
fs.writeFileSync('package.json', outdata);



