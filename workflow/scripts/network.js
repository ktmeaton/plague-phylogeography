var network = {
	nodes:[
		{id:0, country: 'Russia', continent: 'NA', size:89, lat:64.6863136, lon: 97.7453061},
		{id:1, country: 'Lithuania', continent: 'NA', size:4, lat:55.35000030000001, lon: 23.7499997},
		{id:2, country: 'Estonia', continent: 'NA', size:1, lat:58.75237779999999, lon: 25.3319078},
		{id:3, country: 'Germany', continent: 'NA', size:13, lat:51.0834196, lon: 10.4234469},
		{id:4, country: 'China', continent: 'NA', size:233, lat:35.000074, lon: 104.999927},
		{id:5, country: 'Azerbaijan', continent: 'NA', size:16, lat:40.3936294, lon: 47.787250799999995},
		{id:6, country: 'Armenia', continent: 'NA', size:9, lat:40.7696272, lon: 44.6736646},
		{id:7, country: 'Georgia', continent: 'NA', size:9, lat:41.6809707, lon: 44.02873820000001},
		{id:8, country: 'Mongolia', continent: 'NA', size:29, lat:46.825038799999994, lon: 103.8499736},
		{id:9, country: 'Tajikistan', continent: 'NA', size:2, lat:38.62817329999999, lon: 70.8156541},
		{id:10, country: 'Kyrgyzstan', continent: 'NA', size:14, lat:41.5089324, lon: 74.724091},
		{id:11, country: 'England', continent: 'NA', size:9, lat:52.5310214, lon: -1.2649062},
		{id:12, country: 'Spain', continent: 'NA', size:2, lat:39.3260685, lon: -4.8379791},
		{id:13, country: 'France', continent: 'NA', size:8, lat:46.603353999999996, lon: 1.8883334999999999},
		{id:14, country: 'Nepal', continent: 'NA', size:1, lat:28.108392900000002, lon: 84.0917139},
		{id:15, country: 'India', continent: 'NA', size:6, lat:22.3511148, lon: 78.6677428},
		{id:16, country: 'Kazakhstan', continent: 'NA', size:24, lat:47.2286086, lon: 65.20931970000001},
		{id:17, country: 'Turkmenistan', continent: 'NA', size:6, lat:39.37638070000001, lon: 59.392460899999996},
		{id:18, country: 'Iran', continent: 'NA', size:4, lat:32.6475314, lon: 54.56435160000001},
		{id:19, country: 'United States of America', continent: 'NA', size:20, lat:39.7837304, lon: -100.4458825},
		{id:20, country: 'Norway', continent: 'NA', size:1, lat:60.5000209, lon: 9.099971499999999},
		{id:21, country: 'Italy', continent: 'NA', size:1, lat:42.6384261, lon: 12.674297},
		{id:22, country: 'Poland', continent: 'NA', size:1, lat:52.215933, lon: 19.134422},
		{id:23, country: 'Switzerland', continent: 'NA', size:8, lat:46.81333125000001, lon: 8.444947437939406},
		{id:24, country: 'The Netherlands', continent: 'NA', size:2, lat:52.500169799999995, lon: 5.7480820999999995},
		{id:25, country: 'Kenya', continent: 'NA', size:1, lat:-0.16671689999999995, lon: 37.48603689110992},
		{id:26, country: 'Uganda', continent: 'NA', size:1, lat:1.3769751, lon: 32.72391099113675},
		{id:27, country: 'Democratic Republic of the Congo', continent: 'NA', size:2, lat:-2.9814344, lon: 23.82226360000001},
		{id:28, country: 'Indonesia', continent: 'NA', size:3, lat:-2.4833826, lon: 117.89028529999999},
		{id:29, country: 'Madagascar', continent: 'NA', size:3, lat:-18.9249604, lon: 46.441642200000004},
		{id:30, country: 'Vietnam', continent: 'NA', size:5, lat:13.2904027, lon: 108.42651129999999},
		{id:31, country: 'Myanmar', continent: 'NA', size:1, lat:17.1750495, lon: 95.9999652},
		{id:32, country: 'Peru', continent: 'NA', size:65, lat:-6.8699697, lon: -75.0458515},
		{id:33, country: 'Canada', continent: 'NA', size:1, lat:61.0666922, lon: -107.991707},
		{id:34, country: 'Algeria', continent: 'NA', size:1, lat:28.000027199999998, lon: 2.9999825},
		{id:35, country: 'Brazil', continent: 'NA', size:1, lat:-10.333333300000001, lon: -53.2},
		{id:36, country: 'Zimbabwe', continent: 'NA', size:1, lat:-19.01688, lon: 29.353650159713393},
		{id:37, country: 'Bolivia', continent: 'NA', size:1, lat:-17.056869600000002, lon: -64.9912286}
	],
	links:[
		{source:0, target:4, size:5},
		{source:4, target:0, size:5},
		{source:0, target:3, size:3},
		{source:3, target:24, size:1},
		{source:24, target:4, size:1},
		{source:4, target:37, size:1},
		{source:37, target:32, size:1},
		{source:4, target:34, size:1},
		{source:34, target:35, size:1},
		{source:35, target:36, size:1},
		{source:4, target:19, size:1},
		{source:19, target:33, size:1},
		{source:19, target:32, size:1},
		{source:32, target:19, size:1},
		{source:32, target:15, size:1},
		{source:4, target:31, size:1},
		{source:4, target:30, size:1},
		{source:4, target:29, size:1},
		{source:29, target:0, size:1},
		{source:0, target:15, size:1},
		{source:4, target:28, size:1},
		{source:4, target:15, size:3},
		{source:4, target:25, size:1},
		{source:25, target:27, size:1},
		{source:27, target:26, size:1},
		{source:24, target:0, size:1},
		{source:3, target:1, size:1},
		{source:1, target:0, size:1},
		{source:0, target:11, size:1},
		{source:11, target:13, size:1},
		{source:11, target:0, size:1},
		{source:1, target:3, size:2},
		{source:3, target:23, size:1},
		{source:1, target:22, size:1},
		{source:3, target:11, size:2},
		{source:3, target:21, size:1},
		{source:3, target:13, size:2},
		{source:3, target:20, size:1},
		{source:3, target:12, size:2},
		{source:4, target:16, size:1},
		{source:16, target:0, size:13},
		{source:0, target:5, size:3},
		{source:16, target:17, size:6},
		{source:16, target:5, size:4},
		{source:5, target:0, size:1},
		{source:5, target:18, size:1},
		{source:18, target:19, size:1},
		{source:16, target:4, size:1},
		{source:4, target:8, size:4},
		{source:8, target:4, size:3},
		{source:15, target:14, size:1},
		{source:0, target:8, size:2},
		{source:4, target:10, size:5},
		{source:10, target:11, size:1},
		{source:11, target:3, size:1},
		{source:10, target:9, size:1},
		{source:5, target:6, size:3},
		{source:6, target:7, size:2},
		{source:6, target:5, size:1},
		{source:0, target:2, size:1},
		{source:0, target:1, size:1}
	]
}
