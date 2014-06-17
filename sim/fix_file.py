#!/opt/local/bin/python



import sys

import re



for file in sys.argv[1:]:

    fp = open(file)

    L = []

    for line in fp:

        r = re.search("(0\.\d\d\d)\d+(SCOPA|BASIC)", line)

        if r:

            line = re.sub("0\.\d\d\d\d+(SCOPA|BASIC)", r.group(1) + "\t" + r.group(2), line)

        L.append(line)

    fp.close()

    fp = open(file, "w")

    fp.write("\n".join(L))

    fp.close()




