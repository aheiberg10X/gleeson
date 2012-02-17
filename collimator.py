
#A Source is something that the Collimator can operate on
class CollimatorSource :
    #iterator: is something we can call next on, will provide the data
    #          For example, iterating the lines of a file, or rows of a db
    #eqkey: is a function that takes an item from iterator
    #       and returns a sort key. The intent is that heterogenous Sources
    #       will all produce the same kind of key, so Collimator can compare.
    #       ***self.iterator is expected to be sorted with respect to eqkey***
    #integrator: is a function that takes a list of items from iterator
    #            and merges it into some data structure called 'target'
    #allow_absent: Collimator will find the minimum set of items across multipl     
    #              Sources.  Do we ever let this Source be absent from this 
    #              minimum set?
    #group_repeats: if there are multiple entries with the same key, 
    #               do we treat them as one unit and give them all 
    #               to integrator at once?
    def __init__(self, \
                 iterator, \
                 eqkey, \
                 integrator, \
                 allow_absent = True, \
                 group_repeats = True) :
        self.iterator = iterator
        self.eqkey = eqkey
        self.integrator = integrator
        self.allow_absent = allow_absent
        self.group_repeats = group_repeats

# the default way to handle AbsentExceptions
def raiseAE( ae ) : raise ae

#Collimator takes multiple Sources and coalesces their equivalent entries
#Equivalence among the Source entries is defined by self.comparator, which 
#operates on the source.eqkey( source.iteritem ) of two Sources
#comparator - inputs: itemA, itemB
#         output: -1 iff itemA < itemB
#                  0 iff itemA = itemB
#                 +1 iff itemA > itemB
#targetCreator: a function returning a data structure that the Sources expect
#               as the first parameter in their respective integrator functions
class Collimator :
    def __init__(self, \
                 sources, \
                 comparator, \
                 targetCreator, \
                 absentHandler=raiseAE) :
        self.sources = sources
        self.nsources = len(sources)
        self.comparator = comparator
        self.count = 0
        self.targetCreator = targetCreator
        self.absentHandler = absentHandler
        self.exhausted_sources = [False]*self.nsources
        self.next_entry = [self.sources[i].iterator.next() \
                           for i in range(self.nsources)]
        self.entries = [ [] for i in range(self.nsources) ]
        for i in range(self.nsources) :
            self.fillEntry(i)

    #get the next 
    def fillEntry( self, source_ix ) :
        if not self.next_entry[source_ix] :
            self.exhausted_sources[source_ix] = True
            return

        entry = self.next_entry[source_ix]
        key = self.sources[source_ix].eqkey(entry)
        entries = [key, entry]
        try :
            grouping = self.sources[source_ix].group_repeats
            while True :
                next_entry = self.sources[source_ix].iterator.next()
                self.next_entry[source_ix] = next_entry
                if grouping :
                    next_key = self.sources[source_ix].eqkey( next_entry )
                    comp = self.comparator( key, next_key )
                    if comp == 0 :
                        entries.append( next_entry )
                    elif comp == -1 : break
                    else :
                        message = "Source %d has %s before %s" \
                                   % (source_ix, str(key), str(next_key))
                        print message
                        assert False
                else :
                    break

        except StopIteration :
            self.next_entry[source_ix] = False

        self.entries[source_ix] = entries
        #print "source %d new entries" % source_ix, entries

    def __iter__(self) :
        while True :
            non_exhausted = [ i for i in range(self.nsources) \
                                if not self.exhausted_sources[i] ]
            #print "nonexhauted", non_exhausted
            if len(non_exhausted) == 0 :
                raise StopIteration

            current_min = self.entries[non_exhausted[0]][0]
            min_sources = []
            for i in non_exhausted :
                #print "comparing", current_min, "to entrie", i, self.entries[i][0]
                comp = self.comparator( current_min, self.entries[i][0] )
                if comp == 1 :
                    current_min = self.entries[i][0]
                    min_sources = [i]
                elif comp == 0 :
                    min_sources.append(i)

            #print min_sources
            target = self.targetCreator()
            for i in range(self.nsources) :
                if i in min_sources :
                    #print i, self.entries[i][0]
                    target = self.sources[i].integrator( target, \
                                                         self.entries[i][1:] )
                    self.fillEntry(i)
                else :
                    if not self.sources[i].allow_absent :
                        ae = AbsentException(i,current_min,target)
                        self.absentHandler( ae )
            self.count += 1
            yield target

        #except StopIteration :
            #print "All sources have been exhausted"

class AbsentException(Exception) :
    def __init__(self, src_ix, missed_target_key, missed_target) :
        self.ix = src_ix
        self.missed_target_key = missed_target_key
        self.missed_target = missed_target
    def __str__(self) :
        print 'called'
        print self.missed_target_key
        return "Source %d is not allowed to be unmatched with: %s|| " % (self.ix, self.missed_target_key)


def test() : 
    source1 = [1,2,2,3,4]
    source2 = [3,3,3,4,5]
    def eqk( entry ) :
        return entry
    def intg( target, entries ) :
        target.extend(entries)
        return target
    def blankTarget() : return []

    source1 = Source( iter(source1), eqk, intg )
    source2 = Source( iter(source2), eqk, intg )
    c = Collimator( [source1,source2], equiComp, blankTarget )
    for t in c :
        print t
    #c.fillEntry( 0 )
    #print c.entries, c.next_entry
    #c.fillEntry(0)
    #print c.entries, c.next_entry
    #c.fillEntry(1)
    #print c.entries, c.next_entry
##########################################################################
########  Helpers / Simple Defaults   ####################################
##########################################################################

def equiComp( x, y ) : 
    if x < y : return -1
    elif x==y : return 0
    else : return 1

if __name__ == '__main__' :
    test()
