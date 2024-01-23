# Heat simulation of a single heatfin CPU radiator

University assigment in first year of a master's degree in applied mathematics.  
You can find a full report in /doc.


## Compilation

Place yourself in the root of the Project.  
Move or make a fresh build directory and move into it:
```bash
mkdir build
cd build
```
Compile an executable :
```bash
cmake --preset release ..
make
```

## Usage
### Run a simulation

```bash
./run <config_file>
```
Example:
```bash
./run ../simul.cfg
```
**Note:** writing into the files can take some times (with the current parameters there is
around 240'000 * 600 = 144 millions of small lines to write).  

You can modify or make your own **configuration files** as long as it's **exactly** in the 
same **format**.

### Visualization 

For 1d solution you can chose in the `main.cpp` if you want to see the plots or not.    
Additionnaly if you want to see more plots and testings you can use :
```bash
make test
```
Solutions are saved in `data`.  
3d solutions can be visualized in a visualization software like `Paraview`.  
**Note:** to do the plots you will need the `pandas` library which can be installed, using [pip](https://pip.pypa.io/en/stable/)  
```bash
pip install pandas
```
or using [conda](https://www.anaconda.com/)/[mamba](https://github.com/mamba-org/mamba)
```bash
conda install pandas
```  
```bash
mamba install -c conda-forge pandas
```
**Note:** if the data directory gets too big you can delete it, the saves functions knows how to create one without crashing.







