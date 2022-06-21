#!/bin/bash
origin=$1
dest=$2
if [ -f $origin ]; then
   mv $origin $dest
fi
