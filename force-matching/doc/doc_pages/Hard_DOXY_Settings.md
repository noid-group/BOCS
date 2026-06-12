@page Hard_DOXY_Settings Manually Entering Doxygen Settings without a Premade Doxyfile
@author Maria Lesniewski 

This is a page enumerating the settings I used to compile this documentation using the Doxy Wizard GUI for Doxygen. 


The DoxyWizard will open to the main menu, which has 3 main pannels of operation listed below. 
Regenerating the document in the same format as you are viewing it now can be accomplished by linking the appropriate directories 
 and checking the same formatting boxes in each of the 3 main pannels as we describe below. 
The settings that I use can be found inside the Doxyfile file in force-matching/doc/ or pressure-matching/doc/ folders and contains all options listed here, pre loaded.

1)  "WIZARD" Pannel Settings

    1)  PROJECT
        a) Projectname: BOCS
        
        b) Project synopsis: force-matching documentation
        
        c) Project logo: Doctor BOCS log.jpeg

        d) Source code directory: BOCS_Developers-main/force-matching/src
        
        e) [x] Scan Recursively

        f) Destination directory: BOCS_Developers-main/force-matching/doc

        

    2)  MODE
        
        a)    Desired extraction mode: All Entities
        
        b)    Prgramming language to optimize the results for: Optimize for C or PHP output


    3)  OUTPUT
        
        a)     [x] HTML

               (x)   with navigation panel

               [x]   with search function

        *no other boxes checked on this page*

    4)  DIAGRAMS
        
        a)  (o) Use dot tool from the GraphViz package

                [x] Class graphs

                [x] Collaboration diagrams

                [x] Overall Class hierarchy

                [x] Include dependency graphs

                [x] Included by dependency graphs

                [x] Call graphs

                [x] Called by graphs

2)  "EXPERT" Pannel Settings
    
    * If not listed, option is likely left blank or unchecked *

    1)  PROJECT
        a) DOXYFILE_ENCODING: UTF-8
        
        b) PROJECT_NAME: BOCS

        c) PROJECT_NUMBER: v5.0.0

        d) PROJECT_BRIEF: force-matching documentation

        e) PROJECT_LOGO: force-matching/doc/doc_images/Doctor BOCS log.jpeg

        f) OUTPUT_DIRECTORY : force-matching/doc

        g) OUTPUT_LANGUAGE: English

        h) BRIEF_MEMBER_DESC [x]

        i) REPEAT_BRIEF [x]

        j) JAVADOC_AUTOBRIEF [x]

        k) JAVADOC_BANNER [x]

        l) PYTHON_DOCSTRING [x]

        m) INHERIT_DOCS [x]

        n) TAB_SIZE 4

        o) OPTIMIZE_OUTPUT_FOR_C [x]

        p) MARKDOWN_SUPPORT [x]

        q) TOC_INCLUDE_HEADINGS: 5

        r) MARKDOWN_ID_STYLE: DOXYGEN

        s) AUTOLINK_SUPPORT [x]

        t) INLINE_GROUPED_CLASSES [x]

        u) INLINE_SIMPLE_STRUCTS [x]

    2)  BUILD
        
        a) EXTRACT_ALL [x] 

        b) EXTRACT_LOCAL_CLASSES [x] 

        c) RESOLVE_UNNAMED_PARAMS [x]

        d) CASES_SENSE_NAMES SYSTEM

        e) HIDE_SCOPE_NAMES [x]

        f) SHOW_HEADERFILE [x]

        g) SHOW_INCLUDE_FILES [x]

        h) INLINE_INFO [x]

        i) SORT_MEMBER_DOCS [x]

        j) GENERATE_TODOLIST [x]

        k) GENERATE_TESTLIST [x]

        l) GENERATE_BUGLIST [x]

        m) GENERATE_DEPRECATED_LIST

        n) MAX_INITIALIZER_LINES 30

        o) SHOW_FILES [x]

        p) SHOW_NAMESPACES [x]

    3)  MESSAGES
        
        a)  WARNINGS [x]

        b) WARN_IF_UNDOCUMENTED [x]

        c) WARN_IF_DOCERROR [x]

        d) WARN_IF_INCOMPLETE_DOC [x]

        e) WARN_AS_ERROR NO

        f) WARN_FORMAT \$file:\$line:\$text


    4)  INPUT
        a) INPUT: 
                    force-matching/src
                    
                    force-matching/doc/doc_pages
                    
                    force-matching/doc/images_for_doxygen
       
        b) INPUT_ENCODING : UTF-8

        c) RECURSIVE [x]

        d) IMAGE_PATH:
                    force-matching/doc/doc_images
        
        e) FORTRAN_COMMENT_AFTER: 72


    5)  SOURCE BROWSER
        
        a) STRIP_CODE_COMMENT [x]

        b) REFERENCED_BY_RELATION [x]

        c) REFERENCES_LINK_SOURCE [x]

        d) VERBATIM_HEADERS [x]

    6)  INDEX
        
        a) ALPHABETICAL_INDEX [x]  

    7)  HTML
       
        a) GENERATE_HTML [x] 

        b) HTML_OUTPUT: html

        c) HTML_FILE_EXTENSION: .html

        e) HTML_COLOR_STYLE: AUTO_LIGHT
        
        f) HTML_COLORSTYLE_HUE : 220
        
        g) HTML_COLORSTYLE_SAT: 100

        h) HTML_DYNAMIC_MENUS [x]

        i) HTML_CODE_FOLDING [x]

        j) HTML_COPY_CLIPBOARD [x]

        k) HTML_INDEX_NUM_ENTRIES: 100

        l) DISABLE_INDEX [x]

        m) GENERATE_TREEVIEW [x]

        n) FULL_SIDEBAR [x]

        o) ENUM_VALUES_PER_LINE 4

        p) TREEVIEW_WIDTH: 250

        q) OBFUSCATE_EMAILS [x]

        r) HTML_FORMULA_FORMAT png

        s) FORMULA_FONTSIZE: 10

        t) SEARCHENGINE [x]

    8)  LaTeX
   
        * I don't check anything on this page *
  
    9)  RTF
   
        * I don't check anything on this page *

    
    10) Man

         * I don't check anything on this page *
    11) XML
   
    * I don't check anything on this page *
   
    12) DOCBOOK

     * I don't check anything on this page *
    
    13) AUTOGEN

     * I don't check anything on this page *
  
    14) SQLITE3

     * I don't check anything on this page *
  
    15) PERLMOD

     * I don't check anything on this page *

    16) PREPROCESSOR

        a) ENABLE_PREPROCESSING [x]

        b) SEARCH_INCLUDES [x]

        c) SKIP_FUNCTION_MACROS [x]

    17) EXTERNAL
        
        a) EXTERNAL_GROUPS [x]

        b) EXTERNAL_PAGES [x]
  
    18) DOT
   
        a) HIDE_UNDOC_RELATIONS [x]

        b) HAVE_DOT [x]

        c) DOT_COMMON_ATTR : fontname=Helvectica,fontsize=10

        d) DOT_EDGE_ATTR : labelfontname=Helvetica,labelfontsize-10

        e) DOT_NODE_ATTR : shape=box,heigh=0.2,width =0.4

        f) CLASS_GRAPH : YES

        g) COLLABORATION_GRAPH [x]

        h) GROUP_GRAPHS [x]

        i) DOT_WRAP_THRESHOLD 17

        j) INCLUDE_GRAPH [x]

        k) INCLUDED_BY_GRAPH [x]

        l) CALL_GRAPH [x]

        m) CALLER_GRAPH [x]

        n) GRAPHICAL_HIERARCHY [x]

        o) DIRECTORY_GRAPH [x]

        p) DIR_GRAPH_MAX_DEPTH : 1

        q) DOT_IMAGE_FORMAT : .png

        r) DOT_PATH : windows_10_msbuild_Release_graphviz-11.0.0-win32/Graphviz/bin/dot.exe   <-- CHANGE TO YOURS

        s) DOT_GRAPH_MAX_NODES : 50

        t) MAX_DOT_GRAPH_DEPTH : 0

        u) GENERATE_LEGEND [x]
        
        v) DOT_CLEANUP [x]

3)  "RUN" Pannel
   
    No additional checks, just click "Run doxygen"