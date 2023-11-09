# working on testing restart capabilities


# make environment
conda create -n dclaw_mwe python=3.11
conda activate dclaw_mwe
conda install pip

# install a few other dependencies
conda install pyyaml -c conda-forge

# install D-Claw
cd ..
git clone https://github.com/geoflows/D-Claw
export CLAW=.../D-Claw
cd D-Claw/python
pip install -e .

# I ran into issues compiling 
# https://developer.apple.com/forums/thread/739338
# and did the following to fix
# install gfortran from conda-forge
conda install -c conda-forge gfortran
# this messed up my FFLAGS, when I set them back to normal, I compiled without error.
export FFLAGS='-O2 -fopenmp -fdefault-double-8 -fdefault-real-8 -fdefault-integer-8 -fallow-argument-mismatch'

# for posteriority, when compiler flags were messed up, they had been set to 
echo $FFLAGS
-march=core2 -mtune=haswell -ftree-vectorize -fPIC -fstack-protector -O2 -pipe -isystem /opt/anaconda3/envs/dclaw_mwe/include

# run the two step option
cd ../../dam_break_twostep
make new
python setup.py
make .output

# run the restart case
cd ../dam_break_restart
make new
python setup.py
make .data
ln -s ../dam_break_twostep/_output/fort.chk0001 restart.data./xgeoclaw 

# compare the following two files
dam_break_twostep/_output/fort.q0001 <- made without restart
dam_break_restart/fort.q0002 <- made with restart


# additional notes

dam_break_basecase is a basecase of a not so simple mwe for dclaw

to test restart we have two cases 
 
both have only one grid, level 1

- dam_break_twostep, which runs two coarse grid timesteps, generating _output/fort.chk0001. 

- symlink _output/fort.chk0001 to dam_break_restart/restart.data

then run ./xgeoclaw from within that directory, the second coarse timestep will be taken, generating a difference.

We had a hard time getting the make .output to see the restart.data file correctly. 