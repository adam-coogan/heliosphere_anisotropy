#!/usr/bin/env bash

# Use this by running it in the directory with data files to be reformatted.

for runfile in $1/*.csv; do
    firstchar=$(head -c 1 $runfile)

    if [[ $firstchar != '#' ]]; then
        perl -p -i -e 's/E/e/g' $runfile
        perl -p -i -e 's/^\h+//g' $runfile
        perl -p -i -e 's/\h+/,/g' $runfile
        perl -pi -e 'print "# Run exit points.  Columns are r (AU), th (rad), ph (rad), ek (GeV), s (s).  Made with Strauss code. \n" if $. == 1' $runfile
        echo "Reformatted " $runfile
    fi
done

echo "Done reformatting!"


