import globes

a = 0

def memoed() :
    global a
    if a==0 :
        print "unset"
        a = 7
    return a

