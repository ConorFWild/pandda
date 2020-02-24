## Running/testing in Docker container
1. Clone this repository
```
git clone <repository>
```

2. Build the docker container (from the cloned repo directory)
```
docker build -t pandda-2 .
```

3. Run the docker container interactivley
```
docker run -it pandda-2 /bin/bash
```

4. Now inside the container, setup ccp4 env
```
source /ccp4/bin/ccp4.setup-sh
```

5. Test that pandda-2 is set up correctly by using the BAZ2B data
```
ccp4-python /pandda/program/pandda_2_luigi.py data_dirs="/BAZ2B/data/*" pdb_style="*.dimple.pdb" mtz_style="*.dimple.mtz" cpus=8 out_dir="/BAZ2B_out" diffraction_data.structure_factors="2FOFCWT,PH2FOFCWT"
```


## Setup instructions - from source 

1. Create a ccp4 install

2. ```git clone https://github.com/ConorFWild/pandda.git; cd pandda; ccp4-python pip install --upgrade .```


