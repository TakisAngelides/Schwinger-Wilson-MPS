#!/bin/bash

cp ${KRB5CCNAME#FILE:} /tmp/krb5cc_${UID}
unset KRB5CCNAME
condor_submit_dag -f "$@"
exit $?