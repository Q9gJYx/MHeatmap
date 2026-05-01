# data

Store dataset references and lightweight metadata here.

Best practice:

- Avoid committing large/raw/private datasets.
- Keep reproducible acquisition steps in docs.

1. Main synthetic benchmark (Table 1 Synthetic Component) (`main_experiment/synthetic_benchmark`)

- The synthetic benchmark is generated procedurally by code.
- No raw synthetic data files are required or committed.
- Reproducible outputs are written to `output/main_experiment/synthetic_benchmark`.

2. Main real-world benchmark (Teaser, Table 1 Real-world component) (`main_experiment/real_world_benchmark`)

- Raw public data sources are listed under
  `data/main_experiment/real_world_benchmark`.
- Preprocessing and benchmark code are stored under
  `main_experiment/real_world_benchmark`.
- Reproducible outputs are written to `output/main_experiment/real_world_benchmark`.

3. Tabula Sapiens Example (Fig 3) (`examples/tabula_sapiens`)

- This example can be replicated using `examples/tabula_sapiens/tabula_sapien_v2.py` script
- The datasets can be downloaded from: https://figshare.com/articles/dataset/Tabula_Sapien_v2/27921984
- Remember to change the dataset routes in the script to your own paths.

