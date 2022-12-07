#!/usr/bin/env bash
# This script caches the ouput of shell commands when first run and then uses
# the cache if the arguments haven't changed.
# Modified version of https://stackoverflow.com/a/36092050 to avoid race
# condition.

# hash all arguments
KEY="$@"
EXIT_CODE=0

# hash last modified dates of any files
for arg in "$@"
do
  if [ -f "$arg" ]
  then
    KEY+=`date -r "$arg" +\ %s`
  fi
done

# use the hash as a name for temporary file
FILE="/tmp/command_cache.`echo -n "$KEY" | md5sum | cut -c -10`"
FILE_LOCK="${FILE}.lock"

# check file isn't locked
while [ -f $FILE_LOCK ]
do
  sleep 1
done

# use cached file or execute the command and cache it
if [ -f $FILE ]
then
  cat $FILE
  exit $EXIT_CODE
else
  # lock file
  touch $FILE_LOCK
  $@ | tee $FILE
  EXIT_CODE=${PIPESTATUS[0]}
  # delete cached output if command exited with an error
  if [ "$EXIT_CODE" != "0" ]
  then
    rm $FILE
  fi
  rm -f $FILE_LOCK
  exit $EXIT_CODE
fi
