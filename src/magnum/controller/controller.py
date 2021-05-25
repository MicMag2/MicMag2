# Copyright 2012, 2013 by the Micromagnum authors.
#
# This file is part of MicroMagnum.
#
# MicroMagnum is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MicroMagnum is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with MicroMagnum.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import print_function

import multiprocessing as mp

from magnum.logger import logger
from magnum.config import cfg

from .controller_base import ControllerBase

import os

# Controllers defined here: LocalController, PrintParametersController, EnvironmentVariableController, SunGridEngineController

class LocalController(ControllerBase):
    """
    This controller will iterate through all parameter sets, if no --prange
    argument is given. If --prange=a,b is given, the controller will iterate
    through the interval given by range(a,b).
    """

    def __init__(self, run, params, *args, **kwargs):
        super(LocalController, self).__init__(run, params, *args, **kwargs)

        if hasattr(cfg.options, 'prange'):
            prange = cfg.options.prange
        else:
            prange = range(self.num_params)

        self.my_params = []
        for i in prange:
            if i < self.num_params:
                self.my_params.append((i, self.all_params[i]))
            else:
                logger.warn("Ignoring parameter id %s (no such parameter set!)" % i)

        if len(self.my_params) == 0:
            logger.warn("Controller: No parameter sets selected!")

    def start(self):
        if cfg.options.toggle_multi_thread_py or cfg.options.toggle_multi_thread_py_1F:
            threads = []
            output = mp.Queue()
            logger.info("Controller: Using Multithread Call!")
            if cfg.options.toggle_multi_thread_py: logger.info("Controller: Will launch parameter sets in seperate sub-directories!")
            else: logger.info("Controller: Will launch parameter sets in same directory as host script!")

            def helper(idx, param):
                self.logCallMessage(idx, param)
                self.run(*param)
                output.put(1)

            def makeNumStr(num, maxNum, minDigits=3):
                strNum = str(num)
                maxStr = str(maxNum)
                if len(maxStr) < minDigits: lenDiff = 3 - len(strNum)
                else:                       lenDiff = len(maxStr) - len(strNum)
                for addLen in range(lenDiff): strNum = "0" + strNum
                return strNum
    
            for idx, param in self.my_params:
                t = mp.Process(target=helper, args=(idx, param))
                threads.append(t)
            logger.info("Controller: Starting " + str(len(threads)) + " threads.")

            # Start all threads
            if cfg.options.toggle_multi_thread_py:
                cwd = os.getcwd() 
                for x in range(len(threads)):
                    num = makeNumStr(x, len(threads),3)
                    if not os.path.isdir(cwd + "/" + num): os.mkdir(cwd + "/" + num)
                    os.chdir(cwd + "/" + num)
                    threads[x].start()
                    os.chdir(cwd)    
            else:
                for x in threads:
                    x.start()         

            # Wait for all of them to finish
            for x in threads:
                x.join()

        else:
            # start individually in sequence without multi-processing
            for idx, param in self.my_params:
                self.logCallMessage(idx, param)
                self.run(*param)

class PrintParametersController(ControllerBase):
    """
    """

    def __init__(self, run, params, print_num_params, print_all_params, *args, **kwargs):
        super(PrintParametersController, self).__init__(run, params, *args, **kwargs)
        self.print_num_params = print_num_params
        self.print_all_params = print_all_params

    def start(self):
        if self.print_num_params:
            print("NUM_PARAMETERS %s" % self.num_params)
        if self.print_all_params:
            for idx, param in enumerate(self.all_params):
                print("PARAMETER %s %s" % (idx, param))

class EnvironmentVariableController(ControllerBase):
    """
    This controller uses an environment variable to select one parameter set.
    """

    def __init__(self, run, params, env, offset=0, *args, **kwargs):
        super(EnvironmentVariableController, self).__init__(run, params, *args, **kwargs)

        try:
            task_id = int(os.environ[env]) - offset
        except:
            logger.error("Could not read environment variable '%s'." % env)
            raise

        if task_id >= len(self.all_params):
            logger.warn("SGE task id is greater than the number of parameter sets.")
            self.my_params = []
        else:
            self.my_params = [(task_id, self.all_params[task_id])]

    def start(self):
        for idx, param in self.my_params:
            self.logCallMessage(idx, param)
            self.run(*param)

class SunGridEngineController(EnvironmentVariableController):
    """
    This controller uses the 'SGE_TASK_ID' enviroment variable to select one parameter set.
    To be used in conjunction with task arrays using the Sun Grid Engine.
    """

    def __init__(self, run, params, *args, **kwargs):
        super(SunGridEngineController, self).__init__(run, params, env="SGE_TASK_ID", offset=1, *args, **kwargs)
