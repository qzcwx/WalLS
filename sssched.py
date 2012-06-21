#!/usr/bin/env python
import sys
import os
import time
import subprocess
import random

from threading import Thread
from optparse import OptionParser

idleThres = 40

#
# Machines use to run the tasks
#
class Machine:
  def __init__(self, name):
    self.name = name
    self.task = None
    self.isDead = False

  def run(self, task):
    self.task = task
    self.task.machine = self
    self.task.start()

  def reset(self):
    task = self.task
    self.task = None
    return task

  def test(self):
    devnull = open("/dev/null", 'rw')
    pingTest = subprocess.call(["/bin/ping", "-c2", "-w2", "-q", self.name],stdin=devnull, stdout=devnull, stderr=devnull)
    if pingTest == 0:
      sshTest = subprocess.call(["/usr/bin/ssh", "-x", self.name, "true"],stdin=devnull, stdout=devnull, stderr=devnull)
      if sshTest == 0:
        self.isDead = False
        return True
    self.isDead = True
    return False

  def free(self):
    print 'ssh '+self.name+' \"mpstat 1 2 | tail -1 | awk \'{ print \$NF }\' \"'
    p = subprocess.Popen('ssh '+self.name+' \"mpstat 1 2 | tail -1 | awk \'{ print \$NF }\' \"', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
    out = p.stdout.read()
    err = p.stderr.read()
    print 'stderr', err 
    idle =float(out) 
    print self.name,'idle', idle

    if idle > idleThres:
        return True
    else:
        return False

#
# A task is one command executed on one machine
#
class Task(Thread):
  def __init__(self, command):
    Thread.__init__(self)
    self.command = command
    self.machine = None
    self.success = False
    
  def run(self):
    sys.stdout.write(str("%% "+self.machine.name+" running: "+self.command+'\n'))
    sys.stdout.flush()
    sys.stderr.flush()
    exitValue = subprocess.call(["/usr/bin/ssh", "-x", self.machine.name,self.command], stdin=sys.stdin, stdout=sys.stdout, stderr=sys.stderr)
    self.success = (exitValue == 0)

#
# Manage a set of tasks object
#
class Manager(Thread):
  def __init__(self, machinesList, commandsList):
    Thread.__init__(self)

    # List of machines
    machinesTest = {}
    for machine in machinesList:
      machinesTest[machine] = None
      
    self.machinesList = []
    for machine in machinesList:
      lmachine = Machine(machine)

      if machinesTest[machine] is None:
        machinesTest[machine] = lmachine.test()
        if machinesTest[machine] == False:
          sys.stdout.write(str("%% "+lmachine.name+" is apparently dead\n"))
          sys.stdout.flush()
          sys.stderr.flush()
        else:
          sys.stdout.write(str("%% "+lmachine.name+" is ok\n"))
          sys.stdout.flush()
          sys.stderr.flush()

      if machinesTest[machine] == True:  
        self.machinesList.append(lmachine)    
    
    # List of tasks to perform
    self.taskToPerform = []
    for command in commandsList:

      self.taskToPerform.insert(0, Task(command))

  def run(self):
    while len(self.taskToPerform) > 0:
      machine = random.choice(self.machinesList)
      print machine.name, len(self.taskToPerform), 'jobs left'
      if machine.isDead:
          print 'machine is dead'
          machine.test()
      else:
          if machine.free()==True:
            if len(self.taskToPerform) > 0:
              task = self.taskToPerform.pop()
              machine.run(task)
      time.sleep(1.0)
    print 'All jobs submitted'


#
# Fill a set with the content of a file (one item = one non empty line)
#
def readListFile(fileName):
  try:
    f = open(fileName)
  except IOError, (errNo, errMsg):
    sys.exit('I/O error with %s: %s' % (fileName, errMsg))

  result = []
  for line in f:
    strippedLine = line.strip()
    if len(strippedLine) > 0:
      result.append(strippedLine)

  return result



def main():
  # Command line parsing
  cmdLine = OptionParser(description = "Simple script scheduler (sssched) for remote job management")
  cmdLine.add_option("-m", "--machines", 
                      dest = "machinesFileName",
                      action = "store",
                      help = "File with list of machines to use (default:\"machines.lst\")",
                      type = "string")
  cmdLine.add_option("-c", "--commands", 
                      dest = "commandsFileName",
                      action = "store",
                      help = "File with list of commands to perform (default:\"commands.lst\")",
                      type = "string")
  cmdLine.add_option("-b", "--background",
                     dest = "backgroundFileName",
                     action = "store",
                     help = "Run sssched in background, redirecting output into given filename (default: \"\")",
                     type = "string")
  
  (options, args) = cmdLine.parse_args()

  # Get the list of the machines
  if options.machinesFileName is None:
    machinesList = readListFile("machines.lst")
  else:
    machinesList = readListFile(options.machinesFileName)

  # Get the list of the commands
  if options.commandsFileName is None:
    commandsList = readListFile("commands.lst")
  else:
    commandsList = readListFile(options.commandsFileName)

  # Get the list of the commands
  if options.backgroundFileName is not None:
    # do the UNIX double-fork magic, see Stevens' "Advanced
    # Programming in the UNIX Environment" for details (ISBN 0201563177)
    try:
      pid = os.fork()
      if pid > 0:
        # exit first parent
        sys.exit(0)
    except OSError, e:
      print >>sys.stderr, "fork #1 failed: %d (%s)" % (e.errno, e.strerror)
      sys.exit(1)
    # do second fork
    try:
      pid = os.fork()
      if pid > 0:
        # exit from second parent, print eventual PID before
        print "Daemon PID %d" % pid
        sys.exit(0)
    except OSError, e:
      print >>sys.stderr, "fork #2 failed: %d (%s)" % (e.errno, e.strerror)
      sys.exit(1)
    # Connect stdin to /dev/null and stdout/stderr to log file
    devnull = os.open('/dev/null', os.O_RDONLY)
    logfile = os.open(options.backgroundFileName, os.O_CREAT|os.O_RDWR)
    os.dup2(devnull, 0)
    os.dup2(logfile, 1)
    os.dup2(logfile, 2)
    
  # Launch the manager
  manager = Manager(machinesList, commandsList)
  manager.start()
  sys.stdout.flush()
  sys.stderr.flush()

if __name__ == "__main__":
  main()
  sys.exit()
