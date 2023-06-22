#!/usr/bin/env python3

from MAUI_modules import Module_check
from MAUI_modules import Config
from MAUI_modules import PreProcess
from MAUI_modules import Pear
from MAUI_modules import MauiCounting
from MAUI_modules import SeqID2DBid

if __name__ == "__main__":
    Module_check()
    config = Config()
    PreProcess(config)
    Pear(config)
    MauiCounting(config)
    SeqID2DBid(config)
