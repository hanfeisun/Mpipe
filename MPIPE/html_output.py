#! /usr/bin/env python
# Time-stamp: <2011-09-08 15:03:29 sunhf>
"""Module Description: Html template to format the output html

Copyright (c) 2011 Hanfei Sun <hfsun.tju@gmail.com>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  experimental
@version: $Revision$
@author:  Hanfei Sun
@contact: hfsun.tju@gmail.com
"""
def html_tag(tag, cont, *args):
    """
    A html tag maker

    @type  tag: str
    @param tag: tag name
    @type  cont: str
    @param  cont: content between the tag
    @type  args: str
    @param  args: options for the tag

    @rtype:   str
    @return:  Join the strings together to make a tag
    """
    out = "<"+tag
    for key in args:
        out += " "+key[0] +'="'+key[1]+'"'
    out += '>%s</%s>\n'%(cont, tag)
    return out

_make_tag = lambda tag:lambda cont,  *args:html_tag(tag, cont, *args)
_stradd =  lambda *args:reduce(lambda s1, s2:s1+s2, args)

def html_template():
    """
    A html template 

    @rtype:   dict
    @return:  the template with different styles

    """
    
    html = _make_tag("html")
    head = _make_tag("head")
    body = _make_tag("body")
    title = _make_tag("title")
    table = _make_tag("table")
    tr = _make_tag("tr")
    td = _make_tag("td")
    
    page = lambda _title, _content:html(head(title(_title)+body(_content)))
    table_style = [["style", "text-align: left; width: 100%;"], 
                 ["border", "1"], 
                 ["cellpadding", "2"], 
                 ["cellspacing", "2"]]
    first_line_style = [["style", "background-color:#009BC2; font-weight:bold"]]
    second_line_style = [["style", "background-color:#F3F6F7"]]
    
    metric_head_syle = [["style", "background-color:#999999; font-weight:bold"]]
    metric_nice_style_a = [["style", "background-color:#FFFEE8"]]
    metric_nice_style_b = [["style", "background-color:#D9E4F8"]]
    metric_line_syle = [["style", "background-color:#FFFFF"]]
    
    table_ = lambda _content:table(_content, *table_style)
    table_line = lambda col, *styles:(tr("".join(map(td, col)), *styles))
    table_line_first = lambda col:table_line(col, *first_line_style)
    table_line_second = lambda col:table_line(col, *second_line_style)
    table_metric_head = lambda col:table_line(col, *metric_head_syle)
    table_metric_nice_line_a = lambda col:table_line(col, *metric_nice_style_a)
    table_metric_nice_line_b = lambda col:table_line(col, *metric_nice_style_b)
    table_metric_line = lambda col:table_line(col, *metric_line_syle)
    
    return {"page":page, 
            "table":table_, 
            "firstline":table_line_first, 
            "secondline":table_line_second, 
            "line":table_line, 
            "white":table_line, 
            "m_head":table_metric_head,
            "m_nice_line_a":table_metric_nice_line_a,
            "m_nice_line_b":table_metric_nice_line_b,
            "m_line":table_metric_line
            }
