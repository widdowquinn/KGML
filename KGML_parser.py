# (c) The James Hutton Institute 2013
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@hutton.ac.uk
#
# Leighton Pritchard,
# Information and Computing Sciences,
# James Hutton Institute,
# Errol Road,
# Invergowrie,
# Dundee,
# DD6 9LH,
# Scotland,
# UK
# 
# The MIT License
#
# Copyright (c) 2010-2014 The James Hutton Institute
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


""" This module provides classes and functions to parse a KGML pathway map

The KGML pathway map is parsed into the object structure defined in
KGML_Pathway.py in this module.

Classes
KGMLParser             Parses KGML file

Functions
read                   Returns a single Pathway object, using KGMLParser
                       internally
"""

from KGML_pathway import *
try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO
try:
    import xml.etree.cElementTree as ElementTree
except ImportError:
    import xml.etree.ElementTree as ElementTree


def read(handle, debug=0):
    """ Returns a single Pathway object.  There should be one and only
        one pathway in each file, but there may well be pathological
        examples out there.
    """
    iterator = parse(handle, debug)
    try:
        first = iterator.next()
    except StopIteration:
        first = None
    if first is None:
        raise ValueError("No pathways found in handle")
    try:
        second = iterator.next()
    except StopIteration:
        second = None
    if second is not None:
        raise ValueError("More than one pathway found in handle")
    return first


def parse(handle, debug=0):
    """ Returns an iterator over Pathway elements

        handle               file handle to a KGML file for parsing
        debug                integer for amount of debug information
                              to print
        This is a generator for the return of multiple Pathway objects.
    """
    # Check handle
    if not hasattr(handle, 'read'):
        if isinstance(handle, str):
            handle = StringIO(handle)
        else:
            exc_txt = "An XML-containing handle or an XML string " +\
                "must be provided"
            raise Exception(exc_txt)
    # Parse XML and return each Pathway
    for event, elem in \
            ElementTree.iterparse(handle, events=('start', 'end')):
        if event == "end" and elem.tag == "pathway":
            yield KGMLParser(elem).parse()
            elem.clear()


class KGMLParser(object):
    """ Parse a KGML XML Pathway entry into a Pathway object
    """
    def __init__(self, elem):
        self.entry = elem

    def parse(self):
        """ Parse the input elements
        """
        def _parse_pathway(attrib):
            for k, v in attrib.items():
                self.pathway.__setattr__(k, v)

        def _parse_entry(element):
            new_entry = Entry()
            for k, v in element.attrib.items():
                new_entry.__setattr__(k, v)
            for subelement in element.getchildren():
                if subelement.tag == 'graphics':
                    _parse_graphics(subelement, new_entry)
                elif subelement.tag == 'component':
                    _parse_component(subelement, new_entry)
            self.pathway.add_entry(new_entry)

        def _parse_graphics(element, entry):
            new_graphics = Graphics(entry)
            for k, v in element.attrib.items():
                new_graphics.__setattr__(k, v)
            entry.add_graphics(new_graphics)

        def _parse_component(element, entry):
            new_component = Component(entry)
            for k, v in element.attrib.items():
                new_component.__setattr__(k, v)
            entry.add_component(new_component)

        def _parse_reaction(element):
            new_reaction = Reaction()
            for k, v in element.attrib.items():
                new_reaction.__setattr__(k, v)
            for subelement in element.getchildren():
                if subelement.tag == 'substrate':
                    new_reaction.add_substrate(int(subelement.attrib['id']))
                elif subelement.tag == 'product':
                    new_reaction.add_product(int(subelement.attrib['id']))
            self.pathway.add_reaction(new_reaction)

        def _parse_relation(element):
            new_relation = Relation()
            new_relation.entry1 = int(element.attrib['entry1'])
            new_relation.entry2 = int(element.attrib['entry2'])
            new_relation.type = element.attrib['type']
            for subtype in element.getchildren():
                name, value = subtype.attrib['name'], subtype.attrib['value']
                if name in ('compound', 'hidden compound'):
                    new_relation.subtypes.append((name, int(value)))
                else:
                    new_relation.subtypes.append((name, value))
            self.pathway.add_relation(new_relation)

        #==========#
        # Initialise Pathway
        self.pathway = Pathway()
        # Get information about the pathway itself
        _parse_pathway(self.entry.attrib)
        for element in self.entry.getchildren():
            if element.tag == 'entry':
                _parse_entry(element)
            elif element.tag == 'reaction':
                _parse_reaction(element)
            elif element.tag == 'relation':
                _parse_relation(element)
            # Parsing of some elements not implemented - no examples yet
            else:
                # This should warn us of any unimplemented tags
                print "Warning: tag %s not implemented in parser" % element.tag
        return self.pathway


if __name__ == '__main__':
    # Check large metabolism
    pathway = read(open('ko01100.xml', 'rU'))
    print pathway
    for k, v in pathway.entries.items()[:20]:
        print v
    for r in pathway.reactions[:20]:
        print r
    print len(pathway.maps)

    # Check relations
    pathway = read(open('ko_metabolic/ko00010.xml', 'rU'))
    print pathway
    for k, v in pathway.entries.items()[:20]:
        print v
    for r in pathway.reactions[:20]:
        print r
    for r in pathway.relations[:20]:
        print r
    print len(pathway.maps)

    # Check components
    pathway = read(open('ko_metabolic/ko00253.xml', 'rU'))
    print pathway
    for k, v in pathway.entries.items():
        print v
    print len(pathway.maps)

    # Test XML representation
    print pathway.get_KGML()

    # Test bounds of pathway
    print pathway.bounds
