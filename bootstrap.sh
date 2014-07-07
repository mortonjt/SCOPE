#!/bin/sh

aclocal -I m4 --install                    || exit 1
autoheader                                 || exit 1
autoconf                                   || exit 1
automake --add-missing --copy              || exit 1

exit 0
