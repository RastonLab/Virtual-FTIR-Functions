# Virtual-FTIR-Functions

## flask

### app.py

`app.py` contains the Flask endpoints used by the front-end.

### functions.py

`functions.py` holds the functions used by `app.py`. At one time, `functions.py` and `app.py` were one file. However, as functions became longer and more complex. It made more sense to split the two.

## Tests

### cURL

This cURL command has be used in a terminal to test the applications ability to properly return the X and Y coordinates.

```
curl -X POST localhost:5000/post_json \
    -H "Content-type: application/json" \
    -d "{ \
        \"minWave\" : 1900, \
        \"maxWave\" : 2300, \
        \"molecule\" : \"CO\", \
        \"pressure\" : 0.01, \
        \"resolution\" : 1, \
        \"numScan\" : 1, \
        \"zeroFill\" : 0, \
        \"source\" : 3100, \
        \"beamsplitter\" : \"AR_ZnSe\", \
        \"cellWindow\" : \"CaF2\", \
        \"detector\" : \"MCT\" \
    }"
```

### Javascript Fetch

This JavaScript program has been used in a browser console to test the Flask endpoint.

```
function doIt() {
  let params = {
    minWave: 1900,
    maxWave: 2300,
    molecule: "CO",
    pressure: 0.001,
    resolution: 1,
    numScan: 1,
    zeroFill: 0,
    source: 3100,
    beamsplitter: "AR_ZnSe",
    cellWindow: "CaF2",
    detector: "MCT",
  };

  // return fetch("http://ec2-44-203-44-133.compute-1.amazonaws.com/post_json", {
  return fetch("http://localhost:5000/post_json", {
    method: "POST",
    mode: "no-cors",
    cache: "no-cache",
    credentials: "same-origin",
    headers: {
      "Content-Type": "application/json",
    },
    redirect: "follow",
    referrerPolicy: "no-referrer",
    body: JSON.stringify({
      minWave: params.minWave,
      maxWave: params.maxWave,
      molecule: params.molecule,
      pressure: params.pressure,
      resolution: params.resolution,
      numScan: params.numScan,
      zeroFill: params.zeroFill,
      source: params.source,
      beamsplitter: params.beamsplitter,
      cellWindow: params.cellWindow,
      detector: params.detector,
    }),
  });
}

doIt()
  .then((response) => response.text())
  .then((data) => console.log(data));
```
