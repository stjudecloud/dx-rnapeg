#!/bin/bash
# Swaps a symlink and its target.
#
# The procedure is:
# 1. Determine the symlink target file/dir
# 2. rm the symlink
# 3. cp the file (or cp -r the dir) to the old symlink location
# 4. rm the old file/dir
# 5. symlink to the new path from the old path
#
# Note that the symlink is completely resolved to determine the true file
# location.  The consequences of this when there are multiple levels of
# symlinks are illustrated in this example using paths A, B, and C
# Before:
#   A: symlink to B; specified as parameter
#   B: symlink to C
#   C: real file
# After
#   A: real file
#   B: symlink to C (unaltered, but now it's indirect)
#   C: symlink to A
#
# Parameters
# $1 = path to symlink
# $2 = "deep" to follow multiple symlinks

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
ORIGINALLY_SYMLINK=$1

# Validate
if [ ! -e "$ORIGINALLY_SYMLINK" ]
then echo "Symbolic link file does not exist: $ORIGINALLY_SYMLINK" >&2; exit 1
elif [ ! -h "$ORIGINALLY_SYMLINK" ]
then echo "Path is not a symbolic link: $ORIGINALLY_SYMLINK" >&2; exit 1
fi

# Backup symlink
backup=`mktemp`
if ! cp -d $ORIGINALLY_SYMLINK $backup
then echo "Could not backup symlink; no changes were made" >&2; exit 1
fi

# Fully resolve
ORIGINALLY_FILE=`readlink -f "$ORIGINALLY_SYMLINK"`
if [ "$ORIGINALLY_FILE" == "" ]
then echo "Symbolic link target could not be determined: $ORIGINALLY_SYMLINK" >&2; exit 1
elif [ ! -e "$ORIGINALLY_FILE" ]
then echo "Symbolic link target not found: $ORIGINALLY_FILE" >&2; exit 1
elif [ -f "$ORIGINALLY_FILE" ]
then dash_r=
elif [ -d "$ORIGINALLY_FILE" ]
then dash_r="-r"
else echo "Symbolic link target is not a regular file or directory: $ORIGINALLY_FILE" >&2; exit 1
fi

# Do the work
if ! rm $ORIGINALLY_SYMLINK
then echo "Remove of original symlink failed; no changes were made." >&2; exit 1
fi
if ! cp $dash_r $ORIGINALLY_FILE $ORIGINALLY_SYMLINK
then
  if [ -e $ORIGINALLY_SYMLINK ]
  then echo "Copy started but did not complete successfully; target location is unchanged but symlink needs cleanup: $ORIGINALLY_SYMLINK" >&2; exit 1
  else
    if cp -d $backup $ORIGINALLY_SYMLINK
    then echo "Copy failed and symlink was restored; net result is no changes were made" >&2; exit 1
    else echo "Copy failed and could not restore symlink; target location is unchanged but symlink is now missing: $ORIGINALLY_SYMLINK" >&2; exit 1
    fi
  fi
fi
if ! rm $dash_r -f $ORIGINALLY_FILE
then echo "Removal of original file/dir failed; now there are two copies and no symlinks; original location needs cleanup: $ORIGINALLY_FILE" >&2; exit 1
fi
if ! ln -s `readlink -f $ORIGINALLY_SYMLINK` $ORIGINALLY_FILE
then echo "Failed to create symlink in original location; file/dir was copied, but original location is now missing: $ORIGINALLY_FILE" >&2; exit 1
fi

echo Now file: $ORIGINALLY_SYMLINK
echo Now symlink: $ORIGINALLY_FILE
exit 0