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
