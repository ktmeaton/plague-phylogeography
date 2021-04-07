var network = {
    nodes:[
      {i:0, country: "Russia",   continent: "Asia",   size:10, lat:10, lon: 20},
      {i:1, country: "China",    continent: "Asia",   size:81, lat:12, lon: 18},
      {i:2, country: "Mongolia", continent: "Asia",   size:15, lat:11, lon:19},
      {i:3, country: "Germany",  continent: "Europe", size:6, lat:9,  lon: 16},
      {i:4, country: "England",  continent: "Europe", size:4, lat:8,  lon:17},
    ],
    links:[
        {source:1, target:0, size:1},
        {source:2, target:0, size:8},
        {source:0, target:2, size:3},
        {source:3, target:4, size:3},
    ]
}
