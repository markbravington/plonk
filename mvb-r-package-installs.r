# R package installation for CSIRO CKMR course, June 2021

# YOU NEED R 4.0.x or 4.1.x

# ALSO READ & RUN "turn_off_bc.r" probably now, or maybe it's OK to do it after this

# You should install step-by-step as shown, in case something doesn't work (then you'll at least know where it fails).
# *Don't* execute the lines with "stop" on them--- just the ones afterwards.
# If you try to just run the whole script it will crash, deliberately-- because I don't want you to do that.
# Try the installation lines one-at-a-time, and only once you've gotten one to work, ...
# ... should you go on to the next.
# If the installations don't work, it's a basic problem with your R setup (or with CRAN...), which...
# ... you will have to fix yourself I'm afraid

# This first one is only needed because of stupid R version clashes between packages that *I* didn't write :)
install.packages( 'Matrix')

stop( "don't run this script holus bolus; read the damn instructions...")+
install.packages('TMB')

stop( "don't run this script holus bolus; read the damn instructions...")+
install.packages( c( 'mvbutils', 'atease', 'debug', 'offarray', 'deconvodisc'),
    repos=c( 'https://markbravington.github.io/Rmvb-repo', getOption( 'repos')))


# If you've gotten this far without horrible errors then you're proably good to go.
# Nervous people can try the steps below to check further.
# If package 'offarray' works you are probably fine for everything.
# Short test here. Don't bother tryting to understand it; just see if it works
# Output on my system is interspersed below.

# You can copy'n'paste the whole lot below if you like (ie no need to single-step)

library( offarray)

#Warning messages:
#1: replacing previous import ‘utils::help’ by ‘mvbutils::help’ when loading ‘offarray’
#2: replacing previous import ‘utils::?’ by ‘mvbutils::?’ when loading ‘offarray’

test <- offarray( 1:6, first=c( X=3, Y=4), last=c( 5, 5))
test

#      Y
#X      [,4] [,5]
#  [3,]    1    4
#  [4,]    2    5
#  [5,]    3    6

test <- autoloop( A=3:4, B=c( 'yes', 'no'), {A*A})
test

#      B
#A      yes no
#  [3,]   9  9
#  [4,]  16 16

test[4,'no'] <- 0
test

#      B
#A      yes no
#  [3,]   9  9
#  [4,]  16  0
