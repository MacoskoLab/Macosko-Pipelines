#!/bin/bash
set -Eeuo pipefail

gcloud compute instances list --project velina-208320 --filter='status=RUNNING' --format='csv[no-heading](name,networkInterfaces[0].accessConfigs[0].natIP)' > /tmp/gcloud_hosts_1

cat /tmp/gcloud_hosts_* > /tmp/gcloud_hosts
echo "DONE"

# etc hosts modification
cat /etc/hosts | sed -n '/# DELETE AFTER THIS LINE/q;p' | sponge /etc/hosts
echo "# DELETE AFTER THIS LINE FOR GCP! (DON'T WRITE BELOW) " >> /etc/hosts
cat /tmp/gcloud_hosts | awk -F, '{print $2" "$1}' >> /etc/hosts

# ssh config
cat ~/.ssh/config | sed -n '/# DELETE AFTER THIS LINE/q;p' | sponge ~/.ssh/config
echo "# DELETE AFTER THIS LINE FOR GCP! (DON'T WRITE BELOW) " >> ~/.ssh/config
cat /tmp/gcloud_hosts | awk -F, '{print "Host "$1"@    HostName "$2"@    StrictHostKeyChecking no@    User USER@    IdentityFile ~/.ssh/google_compute_engine@"}' | tr '@' '\n' >> ~/.ssh/config

osascript -e 'display notification "done"'
