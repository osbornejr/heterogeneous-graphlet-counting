#!/bin/bash
##script to set up a fresh ec2 instance on RACE hub
export DEBIAN_FRONTEND=noninteractive
#get required packages
sudo apt update && sudo apt install -y dtach kitty mailutils msmtp msmtp-mta curl 
#set permissions for script
#sudo chmod +x get-docker.sh
#install docker using recommended method
#cp studies/swb-s3370543-personal-study/get-docker.sh ./
#sudo ./get-docker.sh
curl -fsSL https://get.docker.com -o get-docker.sh
sudo sh get-docker.sh
#set up rootless docker
sudo sh -eux <<EOF
# Install newuidmap & newgidmap binaries
apt-get install -y uidmap
EOF
dockerd-rootless-setuptool.sh install 
#tidy up docker install
rm get-docker.sh

#uptime email
#setup email to work
cat <<EOF | sudo tee /etc/msmtprc > /dev/null
defaults
auth           on
tls            on
tls_trust_file /etc/ssl/certs/ca-certificates.crt
logfile        /var/log/msmtp.log

account        gmail
host           smtp.gmail.com
port           587
from           rolljobertson@gmail.com
user           rolljobertson@gmail.com
password       brii dlfj hslc kshe
account default : gmail
EOF
sudo chmod 600 /etc/msmtprc
sudo chown root:root /etc/msmtprc
#setup uptime check script
sudo chown -R ec2-user:ec2-user /usr/local/bin
sudo echo '#!/bin/bash
LIMIT=$(( 24 * 60 * 60))
UPTIME=$(awk "{print int(\$1)}" /proc/uptime)
HOST=$(hostname)
if [ "$UPTIME" -ge "$LIMIT" ]; then
  echo "The machine $HOST has been up for more than 24 hours. Current uptime: $UPTIME seconds" | mail -s "Uptime Alert" joelwilliamrobertson@gmail.com
fi' > /usr/local/bin/check_uptime.sh && sudo chmod +x /usr/local/bin/check_uptime.sh
#set cron process to periodically run script
(crontab -l 2>/dev/null; echo "0 * * * * /usr/local/bin/check_uptime.sh") | crontab -

#ENV
#pull repo with submodules
git clone --recurse-submodules https://github.com/osbornejr/heterogeneous-graphlet-counting.git
#change permissions for all dirs in repo
sudo chmod -R 644 heterogeneous-graphlet-counting
sudo chown -R ec2-user:ec2-user heterogeneous-graphlet-counting
#cd into repo
cd heterogeneous-graphlet-counting
#install unison
make unison
