#!/bin/bash

NAME=xtomko02.tar.gz

make pack
mv $NAME ../BACK/KRY2-$(date +%m%d)-$(date +%H%M%S).tar.gz
