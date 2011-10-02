import MySQLdb

endOfIteration = -1

def refineException( message, query ) :
    if "integrity" in message.lower() :
        return SQLDuplicate(message,query)
    else :
        return SQLError(message,query)

class SQLError(Exception) :
    def __init__(self, e, query) :
        self.error = e
        self.query = query
    def __str__(self) :
        return "\nError: %s\n----------------\nQuery: %s\n" \
                % (repr(self.error), self.query)

class SQLDuplicate(SQLError) :
     pass

def catch(func) :
    def inner(*args, **kwargs) :
        try :
            return func(*args, **kwargs)
        except SQLError as mse:
            print mse
            assert False
    return inner

class Conn :
    def __init__(self, switch="localhost") :
        self.connection = self.quickConnect( switch )
        print self.connection
        self.cur = self.connection.cursor()

    def quickConnect( self, switch ) :
        connections = {"gleeson-closet" : \
                           ["132.239.160.134", \
                            "root", \
                            "(umulus88", \
                            "gleeson"], \
                       "localhost" : \
                           ["localhost", \
                            "root", \
                            "pHuc7h35", \
                            "gleeson"] \
                      }
        print connections[switch]
        return self.connect( *connections[switch] )

        #@catch
    def connect( self, host, user, password, db ) :
        return MySQLdb.connect( host=host, \
                                user=user, \
                                passwd=password, \
                                db=db )
    #@catch
    def getColumns(self,table) :
        self.cur.execute("SHOW COLUMNS FROM %s" % table)
        rs = self.cur.fetchall()
        return [r[0] for r in rs]

    #@catch
    def queryScalar( self, query, cast ) :
        try :
            q = self.query(query)
            if len(q) == 0 : return False
            if len(q) > 1 :
                raise SQLError("iqueryScalar getting more than one row")
            w = q[0][0]
            return cast( w )
        except IndexError : 
            return False
        except TypeError : 
            return False

    def getNextID( self, table ) :
        nid = self.queryScalar('''select max(id) from %s''' % table, int)
        if not nid and not nid == 0 : return 0
        else : return nid+1

    #should never return a StopIteration
    def iterate( self, table, cols=['*'], order_by=[] ) :
        columns = ','.join(cols)
        order_clause = ""
        if order_by :
            order_clause = "ORDER BY %s" % ','.join(order_by)
        query =  "select %s from %s %s" % (columns,table,order_clause)
        return self.iterateQuery( query, no_stop = True )

    def iterateQuery(self, query, size=10000, no_stop=False) :
        self.cur.execute( query )

        goon = True
        while goon :
            rows = self.cur.fetchmany( size )
            # print "iSQL: new set of %d rows" % len(rows)
            goon = len(rows) > 0
            if goon :
                for row in rows : yield row
            else :
                if no_stop :
                    print "yielding endOfIteration"
                    while True : yield endOfIteration
                else :
                    raise StopIteration

    #@catch
    def query(self,query) :
        self.cur.execute( query )
        return self.cur.fetchall()

    def queryToFile(self, query, fname ) :
        fout = open(fname,'wb')
        for row in self.query(query) :
            fout.write( '%s\n' % '\t'.join([str(t) for t in row]) )
        fout.close()

    def sanitizeValue( self, value ):
        #quote wrap everything except NULL
        if value == 'NULL' or value == '' : return 'NULL'
        else : return "'%s'" % str(value).replace("'","\\'")

    def insert(self, table, values, columns=[], skip_dupes=False ) :
        values = ','.join( map( self.sanitizeValue, values ) )
        if not columns :
            ins = """INSERT INTO gleeson.%s VALUES( %s );""" \
                  % (table, ','.join(values))
        else :
            columns = ','.join(['`%s`' % c for c in columns])
            ins = "INSERT INTO gleeson.%s (%s) VALUES( %s );" \
                   % (table, columns, values)

        try :
            self.cur.execute( ins )
        except Exception as e :
            exc = refineException( repr(e), ins )
            if not( type(exc) == SQLDuplicate and skip_dupes ) :
                raise exc

    def update( self, table, values, columns, eyeD ) :
        values = map( self.sanitizeValue, values )
        eqpairs = ["%s = %s" % (c,v) for (c,v) in zip(columns,values)]
        string = ', '.join(eqpairs)
        update = '''update %s set %s where id = %d''' % (table,string,eyeD)
        try :
            self.cur.execute( update )
        except Exception as e :
            raise SQLError( e, update )

    def wipe(self, table) :
        self.cur.execute("delete from %s" % table)

    def close(self) :
        self.cursor.close()

if __name__=='__main__' :
    dbc = Conn()
    #dbc.update("Variants",["plateIIII",1],["sourc","chrom"],0)
    r = dbc.query( "select name,description from Genes where id = 0" )
    if not r[0][0] : print "whoopie"

    #c.execute("""
#INSERT INTO  `gleeson`.`Variants` (
#`id` ,
#`chr` ,
#`pos` ,
#`dbSNP` ,
#`ref` ,
#`alt` ,
#`qual` ,
#`filter` ,
#`info` ,
#`source`
#)
#VALUES (
#'78',  'dsf',  '45',  '234',  't',  'h',  '23',  'we',  'erf',  'erf'
#);""")
