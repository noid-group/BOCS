@page   tutorial_page BOCS DOXYGEN TUTORIAL
@subpage Hard_DOXY_Settings
@author Maria Lesniewski

This is a page dedicated to giving Noid lab members the knowledge required to maintain/regenerate this documentation from BOCS src code.
As indicated by the lower banner, this set of .html documents was generated using DOXYGEN, with call graphs from DOT.

# How to Regenerate This Documentation
## Downloading Doxygen and DOT and BOCS 
1) Download Doxygen for your operating system [Doxygen_Home_Page](https://www.doxygen.nl/index.html)

2) DOT is a Graphviz executable supported by Doxygen, Download Graphviz for your operating system [Graphviz_Download_Instructions](https://graphviz.org/download/)

I made this page using the precompiled versions for Windows10/11 and additional markdown files (described more below)

3) Go to the [Noid Group Github](https://github.com/noid-group) and download the current working copy of BOCS including the src/ and doc/ folders.

## Running Doxygen  
To generate this documentation as formatted by me use the following steps:

1) Once downloaded, open DoxyWizard (A front end GUI for Doxygen) as shown below

![](open_doxywizard.png)

2) Locate "Doxyfile" under force-matching (or pressure-matching)/doc. When this file is opened it will have presets analogous to when this documentation was initially constructed. If for some reason Doxygen has changed or you are unable to find "Doxyfile" I have a brute force list of Doxygen settings I used on this page of this tutorial. @ref Hard_DOXY_Settings

![](open_doxyfile.png)

3) Change the Project and Input paths to reflect your directory. 
    In particular, change the highlighted inputs on the WIZARD/PROJECT tab to reflect your desired src/destination directories. 
    [I recommend leaving them relatively the same, but with paths on your PC. 
    If desired, I have also left a copy of the BOCS logo in doc/doc_images]

    Also on the EXPERT/INPUT/ tab alter the "Input" and "Image_Path" inputs analagously. 

![INPUT_PATH_EXAMPLE](doxygen_expert_input_input.png)


![IMAGE_PATH_EXAMPLE](doxygen_expert_input_image_path.png)


Also on the EXPERT/Dot/ tab, make sure "HAVE DOT" is set to yes, and that the path
is pointing to your dot executable in your local install of graphviz. 


![DOT_PATH_EXAMPLE](have_dot.png)    


4) Hit "Run doxygen" on the run tab, check results by navigating to index.html and clicking around resultant docs. 
   

# How to Document the Code to Add More to this Documentation
As we saw above, generating the html format is just a matter of linking source directories and checking the right boxes. 
However, in order to get Doxygen to recognize our documentation inside of the source directory we need to follow certain comment structures. 
## Commenting Functions
Functions in BOCS should be commented according to the following template: 
The /** */ syntax indicates to DOXYGEN that the block of commentary is formatted for doxygen, while the 
@ commands will place the commentary in appropriate portions of doxygen function documentation.
Note the space between the * and the text on each line. 
 
``` 
/**
 * my_integer_returning_funct
 * This is a brief description of my func.
 * This is the more detailed description, it prints array entries of the array I pass via example_array
 * @arg example_array This is a list of doubles I want this function to print
 * @arg len_example This is the length of the array I want this function to print
 * @return The last entry of the example_array passed
 * @note This is an optional note I want to tell people about this function
 * @warn This is an optional warning I want to leave people about this function
 * @see my_integer_returning_funct_helper() // This is a related function I want linked to my function on the docs page
 */
int my_integer_returning_funct(double * example_array, int len_example)
{
 for (int i = 0; i < len_example; i++)
 {
    fprintf(stderr, "Array entry %d: %f \n", i, example_array[i])l
 }
 
 return example_array[len_example - 1] 
}
```
Note that some people prefer to use large banners e.g.
```
/********************************
 * A banner style comment
 ***********************************/
```
This is also good as long as JAVADOC_BANNER is set to yes while compiling documentation.

## Commenting Structures 

Structures should be commented via in line documentation syntax and the \\struct command. 
See the following example

```
/**
 *\struct tW_my_tutorial
 * This is a brief description of my structure
 * This is a more detailed one. 
 * It contains a list of words that describe an excellent example thats part of tW_my_tutorial 
 */
struct tW_my_tutorial
{
    char** an_excellent_example;  /*!< This is how I document in line what an_excellent_example is */
}
```

## Commenting File Heads
So far all we comment are the filename, authors, and a brief description about the file for the files list in the docs.
This commentary appears at the top of .c files. 
e.g. cgff.c begins: 
```
/**
 * @file cgff.c
 * @authors Will Noid, Wayne Mullinax, Joseph Rudzinski, Nicholas Dunn, Michael Delyser
 * @brief MPI driver for the cgff calculation
 *
```

## Adding Markdown Pages 
   Simply adding markdown files headed with the \@page inside of the doc_pages folder will allow you to add your own pages to this documentation. A fast way to learn would be to click the doc pages (e.g. open force-matching/doc/doc_pages/tutorial_page.md) and compare to the output .html syntax you are reading. 

   For more information on Markdown: [Markdown_Guide] (https://www.markdownguide.org/getting-started/#additional-resources)
   For a tutorial on Doxygen / Markdown Page Linking: [Doxygen_Basics on Youtube](https://www.youtube.com/watch?v=TtRn3HsOm1s)