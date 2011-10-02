import globes
import time
class Annotator :
    #'name' must match one of the dirs in globes.INT_DIR
    def __init__( self, name ) :
        self.root_dir = "%s/%s" % (globes.INT_DIR, name)
        self.name = name
        self.iterator = False
        self.comp = False,
        self.ltfunc = False,
        self.eqfunc = False,
        self.gtfunc = False,
        self.fast_forward_one_to_many = False
        self.allow_unmatched = False


    #takes a varList (probably built from a Broad VCF file) and writes a file that will
    #be used as input to the actual annotator software
    def parse(self, varList) : pass

    #take the input file made by parse and run it through the annotator
    #output the annotations
    def run(self, args) : pass

    # Preparing for merge
    # Set: self.iterator (turn the output file from run() into an iterator)
    #      self.comp, ltfunc, eqfunc, gtfunc work on (a,b), where 'a' is a member
    #      of iTargets and 'b' is a member of self.iterator
    # dargs is a dictionary of args passed in to let register() know how set these things 
    def register(self, dargs) : pass

    #take the annotations from in self.iterator() and merge them into the targets in iTarget
    def annotate(self, iTargets) :
        w = Walker( iTargets, \
                    self.iterator, \
                    self.comp, \
                    self.ltfunc, \
                    self.eqfunc, \
                    self.gtfunc, \
                    self.fast_forward_one_to_many, \
                    self.allow_iter1_unmatched, \
                    self.allow_iter2_unmatched )
        w.walk()


def equiComp( x, y ) : 
    if x < y : return -1
    elif x==y : return 0
    else : return 1

def announce( x, y ) : print "%s == %s" % (x,y)
def doNothing( x, y ) : pass

# Walker processes two iteratables that are sorted with respect to comp()
# comp() :  Input : item from iter1 'a', items from iter2 'b'
#           Return : -1 if 'a' < 'b', 0 if 'a'=='b', 1 if 'a' > 'b'
# The xxfunc(a,b) functions are the actions to be taken for each of comp()'s 
# return values. They do not return anything, rather they have a side-effect
# that updates one or both of the input items.
# walker() supports one-to-one and one-to-many relationships between iter1 and iter2
class Unmatched1(Exception) : pass
class Unmatched2(Exception) : pass
class Unmatched(Exception) : pass

class AnnoStopIteration( StopIteration ) : pass
class TargetStopIteration( StopIteration ) : pass

def log( func ) :
    def inner( *args, **kwargs ) :
        thing = func( *args, **kwargs )
        #args[0].fout.write( "\n|%s|\n" % thing )
        if thing : args[0].fout.write( "%s\n" % str(thing) )
    return inner

#A walker with one target and mulitple annotators
#once a target item is annotated it is passed to a callback
#Relies on objects coming in to be have certain items (see Annotator)
#annos is [Annotatable, Annotator, Annotator, ... ]
class Walker2 :
    def __init__(self, annos, name="TR") :
        self.annos = annos
        for a in annos : print type(a)
        self.nannos = len(annos)
        #target keeps track of all annotators who have matched it
        self.unmatched = [[False]*(self.nannos-1)] + [False]*(self.nannos-1)
        self.exhausted = [False]*self.nannos
        self.values = [-1]*len(annos)
        self.eq_counts = [0]*len(annos)
        self.fout = open("%s_messages.txt" % name, 'wb')

    def __del__(self) : self.fout.close()

    @log
    def gtfunc( self, annoix, a, b ) :
        #do something?
        return self.annos[annoix].gtfunc(a,b)

    @log
    def eqfunc( self, annoix, a, b ) :
        self.eq_counts[annoix] += 1
        return self.annos[annoix].eqfunc(a,b)

    @log
    def ltfunc( self, annoix, a, b ) :
        #do something?
        return self.annos[annoix].ltfunc(a,b)

    @log
    def getNext(self,ix) :
        message = ""
        if ix == 0 :
            num_unmatched = sum( self.unmatched[0] )
            if num_unmatched != 0 and not self.annos[ix].allow_unmatched :
                message = "annoix: %d, value |%s|, unmatched annos: %s" \
                           % (ix,self.values[ix], str(self.unmatched[0]) )
                print message
                assert False
            self.unmatched[ix] = [True]*(self.nannos-1)
            try :
                self.values[ix] = self.annos[ix].iterator.next()
                return message
            except StopIteration :
                self.exhausted[ix] = True
                raise TargetStopIteration()
        else :
            if self.unmatched[ix] and not self.annos[ix].allow_unmatched :
                message = "annoix: %d, value: |%s| unmatched" \
                             % (ix,self.values[ix]) 
                print message
                assert False
            self.unmatched[ix] = True
            try :
                self.values[ix] = self.annos[ix].iterator.next()
                return message
            except StopIteration :
                self.exhausted[ix] = True
                raise AnnoStopIteration()

    def iterate(self) :
        nxts = [self.getNext(i) for i in range(self.nannos)]
        while True :
            target = self.values[0]
            non_exhausted = [ix for ix in range(1,self.nannos) \
                                       if not self.exhausted[ix] ]
            for annoix in non_exhausted :
                try :
                    anno = self.annos[annoix]

                    #walk each annotator to equality or greater
                    while True :
                        eq = anno.comp( target, self.values[annoix] )
                        if not eq == 1 : break
                        else :
                            self.gtfunc( annoix, target, self.values[annoix] )
                            self.getNext( annoix )

                    #if/when it finally matches the target,
                    #handle the one-to-many case
                    if eq == 0 :
                        self.unmatched[annoix] = False
                        self.unmatched[0][annoix-1] = False

                        self.eqfunc( annoix, target, self.values[annoix] )
                        #one to many here!
                        self.getNext(annoix)
                        nexteq = anno.comp( target, self.values[annoix] )
                        while nexteq == 0 :
                            self.unmatched[annoix] = False
                            self.eqfunc( annoix, target, self.values[annoix] )
                            self.getNext(annoix)
                            nexteq = anno.comp( target, self.values[annoix] )
                        assert nexteq == -1

                    elif eq == -1 :
                        self.ltfunc( annoix, target, self.values[annoix] )
                except AnnoStopIteration :
                    continue

            #apply the callback to the target
            yield self.values[0]
            #move on
            try :
                self.getNext(0)
            except TargetStopIteration :
                raise StopIteration()

#Walker class that checked for skipped items
class Walker :
    def __init__(self, iter1, iter2, \
                 comp =   equiComp, \
                 ltfunc = doNothing, \
                 eqfunc = doNothing, \
                 gtfunc = doNothing, \
                 fast_forward_one_to_many = False, \
                 allow_iter1_unmatched = True, \
                 allow_iter2_unmatched = True) :
        self.iter1 = iter1
        self.iter2 = iter2
        self.comp = comp
        self.ltfunc = ltfunc
        self.eqfunc = eqfunc
        self.gtfunc = gtfunc
        self.fast_forward_one_to_many = fast_forward_one_to_many
        self.allow_iter1_unmatched = allow_iter1_unmatched
        self.allow_iter2_unmatched = allow_iter2_unmatched

        self.iter1_unmatched = False
        self.iter2_unmatched = False

    def getNext1(self,old) :
        if self.iter1_unmatched and not self.allow_iter1_unmatched :
            print "unmatched t1: %s" % old
            #raise Unmatched1()
        self.iter1_unmatched = True
        return self.iter1.next()

    def getNext2(self,old) :
        if self.iter2_unmatched and not self.allow_iter2_unmatched :
            print "unmatched t2 %s" % old
            #raise Unmatched2()
        self.iter2_unmatched = True
        return self.iter2.next()

    def walk(self) :
        t1,t2 = self.getNext1(""), self.getNext2("")
        matches = 0
        try :
            while True :
                eqtest = self.comp(t1,t2)
                if eqtest == -1 :
                    self.ltfunc(t1,t2)
                    t1 = self.getNext1(t1)
                elif eqtest == 0 :
                    self.eqfunc( t1,t2 )
                    matches += 1

                    self.iter1_unmatched = False
                    self.iter2_unmatched = False
                    #handle one-to-many relationships
                    t2 = self.getNext2(t2)
                    #option to fast-forward, not reapplying eqfunc for each t2 in the 
                    #matching run
                    if self.fast_forward_one_to_many :
                        while self.comp(t1,t2) == 0 :
                            t2 = self.getNext2(t2)
                        t1 = self.getNext1(t1)
                        next_eqtest = self.comp(t1,t2)
                    else :
                        next_eqtest = self.comp(t1,t2)
                        if next_eqtest == -1 : t1 = self.getNext1(t1)
                        elif next_eqtest == 0 : 
                            #print "one to many: t2 - %s" % t2
                            #assert False
                            pass
                    if next_eqtest == 1:
                        print t1, t2.getPosition()
                        raise Exception("error in ordering")

                else :
                    self.gtfunc( t1,t2 )
                    t2 = self.getNext2(t2)

        except Unmatched1 :
            s="Unmatched t1: %s" % t1
            raise Exception(s)
        except Unmatched2 :
            s = "Unmatched t2: %s" % t2
            raise Exception(s)
        except StopIteration :
            print "%d matches found" % matches

#first walker, take two iterables and match em up
def walker( iter1, iter2, \
            comp =   equiComp, \
            ltfunc = doNothing, \
            eqfunc = doNothing, \
            gtfunc = doNothing, \
            fast_forward_one_to_many = False) :
    t1,t2 = iter1.next(), iter2.next()
    matches = 0
    try :
        while True :    
            eqtest = comp(t1,t2)
            if eqtest == -1 :
                ltfunc(t1,t2)
                t1 = iter1.next()
            elif eqtest == 0 :
                matches += 1
                eqfunc( t1,t2 )

                #handle one-to-many relationships
                t2 = iter2.next()
                #option to fast-forward, not reapplying eqfunc for each t2 in the 
                #matching run
                if fast_forward_one_to_many :
                    while comp(t1,t2) == 0 :
                        t2 = iter2.next()
                    t1 = iter1.next()
                    next_eqtest = comp(t1,t2)
                else :
                    next_eqtest = comp(t1,t2)
                    if next_eqtest == -1 : t1 = iter1.next()
                    elif next_eqtest == 0 : pass

                if next_eqtest == 1:
                    print t1, t2.getPosition()
                    raise Exception("error in ordering")

            else :
                gtfunc( t1,t2 )
                t2 = iter2.next()
    except StopIteration :
        print "%d matches found" % matches


def gt(a,b) : print "gt: %s - %s" % (a,b)
def lt(a,b) : print "lt: %s - %s" % (a,b)
def eq(a,b) : print "eq: %s - %s" % (a,b)

class WalkerTester :
    def __init__(self,it) :
        self.iterator = it
        self.gtfunc = gt
        self.ltfunc = lt
        self.eqfunc = eq
        self.comp = globes.compareHelper
        self.allow_unmatched = False


if __name__=='__main__' :
    o = iter([1,2,4,8,9])
    t = iter([1,1,2,2,2,4,7,9,10])
    h = iter([4,8])
    def callback( t ) :
        print "||%s||\n" % t

    target = WalkerTester(o)
    anno1 = WalkerTester(t)
    anno2 = WalkerTester(h)
    w = Walker2( [target,anno1,anno2], "Walker2Test" )
    for t in w.iterate() :
        callback(t)

    #w = Walker( o,t,equiComp,doNothing,announce,doNothing,
                #False,False,True )
    #w.walk()



