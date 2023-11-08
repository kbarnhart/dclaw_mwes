# working on testing restart capabilities


# make environment
conda create -n dclaw_mwe
conda activate dclaw_mwe
conda install pip

# install D-Claw
git clone https://github.com/geoflows/D-Claw
export CLAW=.../D-Claw
cd D-Claw/python
pip install -e .

# run the two step option
cd ../../dam_break_twostep
make new
make .output

# run the restart case
cd ../dam_break_restart
make new
make .data
ln -s ../dam_break_twostep/_output/fort.chk0001 restart.data
./xgeoclaw 

# compare the following two files
dam_break_twostep/_output/fort.q0001 <- made without restart
dam_break_restart/fort.q0002 <- made with restart


dam_break_basecase is a basecase of a not so simple mwe for dclaw

to test restart we have two cases 
 
both have only one grid, level 1

- dam_break_twostep, which runs two coarse grid timesteps, generating _output/fort.chk0001. 

- symlink _output/fort.chk0001 to dam_break_restart/restart.data

then run ./xgeoclaw from within that directory, the second coarse timestep will be taken, generating a difference.

We had a hard time getting the make .output to see the restart.data file correctly. 