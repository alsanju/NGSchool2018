# NGSchool2018

## Starting point

Login to remote server:

```
ssh username@147.228.242.14
```

You would need to configure conda environments. If you have already dowloaded miniconda in a previous workshop, just type:

```
bash
source activate albasj
```

If not, you would need to run the following commands:

```
# download miniconda
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
# install miniconda
bash Miniconda2-latest-Linux-x86_64.sh
# add path to the pre-configured shared environments
# create config file in the home folder
touch .condarc
# put the following text into the file and save (use your favorite text editor)
envs_dirs:
 - /mnt/shared_conda/envs
# re-bash
bash
# check if everything is ok and you see the envs
conda info -e
# activate needed env
source activate albasj
```

Now, go to the workshop directory:

```
cd /mnt/albasj/analysis
```

Create your own directory, and change to the dir:

```
mkdir username
cd username
```
