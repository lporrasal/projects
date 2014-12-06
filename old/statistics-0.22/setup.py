#!/usr/bin/env python

import sys, os, time
from distutils.core import setup, Extension
from distutils.command.config import config
try:
  from distutils.command.config import log
except:
  pass

import numpy

run_config = 0
for keyword in sys.argv:
    if keyword=='config':
        try:
            if sys.version_info[:2] < (2, 3):
                print('WARNING: If the configuration fails, please update to Python 2.3 or later')
        except AttributeError:
            print('ERROR: Your version of Python is very old. Please update to Python 2.3 or later')
            sys.exit()
        run_config = 1
        path = os.environ.get("PATH")
        if path:
            directories = path.split(':')
        else:
            directories = []
        directories.insert(0, '.')
        path = ':'.join(directories)
        os.environ['PATH'] = path


extra_compile_args = []
extra_link_args = []
if sys.platform != 'darwin':
    extra_link_args = ['-s']

class config_math (config):

    def run (self):
        try: log.set_verbosity(0)
        except: pass
        self.dump_source = 0
        filename = os.path.join("src", "config.h")
        self.configfile = open(filename,'w')
        self.configfile.write('/* config file from setup.py script ' + time.ctime() + '\n')
        self.configfile.write(' * run on ' + sys.platform + ' */\n')
        self.config_toplevel()
        self.configfile.close()
        print('wrote config file (%s)' % filename)

#----------------------------------------------------------------------
    def config_toplevel(self):
        print("  ============= begin statistics configuration ================")
        # check erf presence, emulate otherwise
        testcode = """\
#include <math.h> 
int main(void)
{ double x=erf(1.);
  if (x<0.842700 || x>0.842702) return 1;
  else return 0;
}
"""
        if self.try_run(testcode,libraries=["m"]):
            print("using erf found in libm")
            self.configfile.write("#define HAVE_ERF\n")
        else:
            print("libm does not contain erf; skipping this feature")

        # check erfc presence, emulate otherwise
        testcode = """\
#include <math.h> 
int main(void)
{ double x=erfc(1.);
  if (x<0.157299 || x>0.157300) return 1;
  else return 0;
}
"""
        if self.try_run(testcode,libraries=["m"]):
            print("using erfc found in libm")
            self.configfile.write("#define HAVE_ERFC\n")
        else:
            print("libm does not contain erfc; skipping this feature")

        print("  ============= end statistics configuration ==================")

if not run_config:
    filename = os.path.join("src", "config.h")
    try:
        file = open(filename)
        for line in file:
            if line.startswith("#define"):
                words = line.split()
                argument = "-D" + words[1]
                extra_compile_args.append(argument)
        file.close()
    except:
        pass
	


extension = Extension("statistics",
                      ["src/statistics.c", "src/statisticsmodule.c"],
                      include_dirs=['src',numpy.get_include()],
                      extra_compile_args=extra_compile_args,
                      extra_link_args=extra_link_args)


setup(name="statistics",
      version="0.22",
      description="Statistics for Python",
      author="Michiel de Hoon",
      author_email="mdehoon@gsc.riken.jp",
      url="http://bonsai.ims.u-tokyo.ac.jp/~mdehoon/software/software.html",
      license="Python License",
      ext_modules=[extension],
      cmdclass = {'config': config_math}
      )
