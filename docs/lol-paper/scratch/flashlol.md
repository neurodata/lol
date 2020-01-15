# STEPS

1. AMI Type

`ami-0698bcaf8bd9ef56d`

Connect:

```
ssh -i /path/to/pemkey ubuntu@<ip>
```

2. Download FlashLOL docker:

```
docker pull ericw95/flashlol:0.0.1
```

3. Format and Drives:

```
# format drives
sudo mkfs.ext3 /dev/xvdb
sudo mkfs.ext3 /dev/xvdc

# mount drives
sudo mkdir -p /mnt/ssd1
sudo mkdir -p /mnt/ssd2

sudo mount /dev/xvdb /mnt/ssd1
sudo mount /dev/xvdc /mnt/ssd2
```

4. Download SWU4 data and make available:

```
screen -R
sudo aws s3 cp --no-sign-request --recursive s3://mrneurodata/data/SWU4/ndmg_0-0-48/reg_dti/ /mnt/ssd1
chmod -R 777 /mnt/ssd1
```

5. Clone repo:

```
mkdir ~/Documents
cd ~/Documents
git clone https://github.com/neurodata/lol.git
```

6. Launch Docker:

```
docker run -ti --entrypoint /bin/bash -v /home/ubuntu/Documents/lol:/lol -v /mnt/ssd1:/brains ericw95/flashlol:0.0.1
```

7. Start memory logging:

```
cd ~/Documents/lol/docs/lol-paper/scratch
nohup ./memlog.sh > mem.txt &
```

8. Run scripts:

```
cd /lol/docs/lol-paper/scratch/
```
