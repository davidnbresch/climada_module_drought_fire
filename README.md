# climada_module_drought_fire
climada module for drought and (bush)fire hazard

Purpose

This module implememts drought and (bush)fire hazards. 

Get to know climada

    Go to the wiki and read the introduction and find out what climada and ECA is.
    Are you ready to start adapting? This wiki page helps you to get started!
    Read more on natural catastrophe modelling and look at the GUI that we have prepared for you.
    Read the core climada manual (PDF).

Set-up

In order to grant core climada access to additional modules, create a folder ‘climada_modules’ on the same level your core climada folder resides and copy/move any additional modules into climada_modules. You might shorten the module filename(s), i.e. without 'drought_fire' instead of 'climada_module_drought_fire' as a folder name.

E.g. if the addition module is named climada_module_MODULE_NAME, we should have

    .../climada the core climada, with sub-folders as
    .../climada/code
    .../climada/data
    .../climada/docs

and then

    .../climada_modules/MODULE_NAME with sub-folders as
    .../climada_modules/MODULE_NAME/code
    .../climada_modules/MODULE_NAME/data
    .../climada_modules/MODULE_NAME/docs

this way, climada sources all modules' code upon startup

copyright (c) 2016, David N. Bresch, david.bresch@gmail.com all rights reserved.
