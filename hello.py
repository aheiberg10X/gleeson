#!/usr/bin/python
import md5
import cgi
import explanations
import sys
sys.path.append('/home/Gleeson/variant_filtering/andrew')
print "Content-Type: text/html"     # HTML is following
print                               # blank line, end of headers
print sys.path
#import query


#fields = cgi.FieldStorage()
#if 'passwd' in fields :
    #m = md5.new()
    #m.update(fields['passwd'].value)
    #hashed = m.hexdigest()
    #if hashed == "3712907b986216d86a7598e1912cc405" :
        #print '<form name="query" action="hello.py">'
        #print '    <p><select name="interval">'
        #print '        <option value="Gene">Gene</option>'
        #print '        <option value="Exon">Exon</option>'
        #print '        <option value="Variant">Variant</option>'
        #print '    </select></p>'
        #print '    <p><input type="text" name="include_atleast" />must be an integer</p>'
        #print '    <p><select name="include_genotype" />'
        #print '        <option value="Hom">Hom</option>'
        #print '        <option value="Het">Het</option>'
        #print '        <option value="Both">Both</option>'
        #print '    </select></p>'
        #print '    <p><input type="text" name="group" />must be a comma separated list of patient_ids, i.e ''200,202,203''. See below on how to lookup these ID''s</p>'
        #print '    <p><select name="exclude_genotype" />'
        #print '        <option value="Hom">Hom</option>'
        #print '        <option value="Het">Het</option>'
        #print '        <option value="Both">Both</option>'
        #print '    </select></p>'
        #print '    <p><input type="text" name="exclude_group" />must be a comma separated list of patient_ids, i.e ''201,204''. See below on how to lookup these ID''s</p>'
        #print '    <p><input type="text" name="call_depth" />must be an integer.  Can be empty</p>'
        #print '    <p><input type="text" name="call_qual" />must be an integer.  Can be empty</p>'
        #print '    <p><input type="text" name="variant_restrictions" />edit this at your own risk</p>'
        #print '    <p><input type="submit" value="Run Query" /></p>'
        #print '</form>'

        #print explanations.template
    #else :
        #print '<p>get lost</p>'
#elif 'interval' in fields :
    ##make params dict
    #params = {}
    #print query.makeReport( params ) 

#else :
    #print "nothing to see here"




