Helpful links
---------------
* [Workflow Description Language](https://github.com/openwdl/wdl/blob/legacy/versions/1.0/SPEC.md)
* [fiss/firecloud/api.py](https://github.com/broadinstitute/fiss/blob/master/firecloud/api.py)

gcloud commands
---------------
```
gcloud auth application-default login --scopes 'openid,https://www.googleapis.com/auth/userinfo.email,https://www.googleapis.com/auth/cloud-platform,https://www.googleapis.com/auth/sqlservice.login,https://www.googleapis.com/auth/drive,https://www.googleapis.com/auth/spreadsheets'
gcloud auth application-default set-quota-project velina-208320

```
```
gcloud auth login
gcloud config set project velina-208320
```

.h5 commands
------------
* `h5dump -n SBcounts.h5`: list all contents
* `h5dump -d /metadata/num_reads SBcounts.h5`: print a specific dataset
