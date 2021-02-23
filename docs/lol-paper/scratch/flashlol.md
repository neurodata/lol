# STEPS

1. AMI Type

`ami-0698bcaf8bd9ef56d`

Connect:

```
ssh -i /path/to/pemkey ubuntu@<ip>
```

2. Download FlashLOL docker:

```
Dataset="SWU4"
docker pull ericw95/flashlol:0.0.3
```

3. Format and Drives:

```
# format drives
sudo mkfs.ext3 /dev/xvdb

# mount drives
sudo mkdir -p /mnt/ssd1

sudo mount /dev/xvdb /mnt/ssd1

sudo mkdir -p /mnt/ssd1/dwi
```

4. Download Data data and make available:

```
screen -R
sudo aws s3 cp --no-sign-request --recursive s3://mrneurodata/data/$Dataset/ndmg_0-0-48/reg_dti/ /mnt/ssd1/dwi
sudo aws s3 cp --no-sign-request --recursive s3://mrneurodata/data/$Dataset/ndmg_0-0-48/reg_dti/ /mnt/ssd1/
sudo chmod -R 777 /mnt/ssd1
```

5. Clone repo:

```
mkdir ~/Documents
cd ~/Documents
git clone https://github.com/neurodata/lol.git
```

6. Launch Docker:

```
docker run -ti --entrypoint /bin/bash -v /<path>/<to>/<repo>/lol:/lol -v /mnt/ssd1:/brains ericw95/flashlol:0.0.1
```

8. Run scripts:

```
cd /lol/docs/lol-paper/scratch/
```
