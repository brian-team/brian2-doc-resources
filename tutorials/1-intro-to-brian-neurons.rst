Introduction to Brian part 1: Neurons
=====================================


.. only:: html

    .. |launchbinder| image:: http://mybinder.org/badge.svg
    .. _launchbinder: https://mybinder.org/v2/gh/brian-team/brian2-binder/master?filepath=tutorials/1-intro-to-brian-neurons.ipynb

    .. note::
       This tutorial is a static non-editable version. You can launch an
       interactive, editable version without installing any local files
       using the Binder service (although note that at some times this
       may be slow or fail to open): |launchbinder|_

       Alternatively, you can download a copy of the notebook file
       to use locally: :download:`1-intro-to-brian-neurons.ipynb`

       See the :doc:`tutorial overview page <index>` for more details.



All Brian scripts start with the following. If you’re trying this
notebook out in the Jupyter notebook, you should start by running this
cell.

.. code:: ipython2

    from brian2 import *


.. parsed-literal::

    /mnt/data/anaconda2/lib/python2.7/site-packages/pkg_resources/py2_warn.py:19: UserWarning: ************************************************************
    You are running Setuptools on Python 2, which is no longer
    supported and
    >>> SETUPTOOLS WILL STOP WORKING <<<
    in a subsequent release. Please ensure you are installing
    Setuptools using pip 9.x or later or pin to `setuptools<45`
    in your environment.
    If you have done those things and are still encountering
    this message, please comment in
    https://github.com/pypa/setuptools/issues/1458
    about the steps that led to this unsupported combination.
    ************************************************************
      sys.version_info < (3,) and warnings.warn("*" * 60 + msg + "*" * 60)


::


    

    ImportErrorTraceback (most recent call last)

    <ipython-input-1-0662cb4848fc> in <module>()
    ----> 1 from brian2 import *
    

    /home/marcel/programming/brian2/brian2/__init__.py in <module>()
         64 __release_date__ = '2019-12-20'
         65 
    ---> 66 from brian2.only import *
         67 from brian2.only import test
         68 


    /home/marcel/programming/brian2/brian2/only.py in <module>()
         13 
         14 from brian2.units import *
    ---> 15 from brian2.utils import *
         16 from brian2.core.tracking import *
         17 from brian2.core.names import *


    /home/marcel/programming/brian2/brian2/utils/__init__.py in <module>()
          3 '''
          4 
    ----> 5 from .logger import *
          6 
          7 __all__ = ['get_logger', 'BrianLogger', 'std_silent']


    /home/marcel/programming/brian2/brian2/utils/logger.py in <module>()
         24 
         25 import brian2
    ---> 26 from brian2.core.preferences import prefs, BrianPreference
         27 
         28 from .environment import running_from_ipython


    /home/marcel/programming/brian2/brian2/core/preferences.py in <module>()
          6 import re
          7 import os
    ----> 8 from collections.abc import MutableMapping
          9 from io import BytesIO
         10 


    ImportError: No module named abc


Later we’ll do some plotting in the notebook, so we activate inline
plotting in the notebook by doing this:

.. code:: ipython2

    %matplotlib inline

If you are not using the Jupyter notebook to run this example (e.g. you
are using a standard Python terminal, or you copy&paste these example
into an editor and run them as a script), then plots will not
automatically be displayed. In this case, call the ``show()`` command
explicitly after the plotting commands.

Units system
------------

Brian has a system for using quantities with physical dimensions:

.. code:: ipython2

    20*volt


::


    

    NameErrorTraceback (most recent call last)

    <ipython-input-3-cd04fc03abe9> in <module>()
    ----> 1 20*volt
    

    NameError: name 'volt' is not defined


All of the basic SI units can be used (volt, amp, etc.) along with all
the standard prefixes (m=milli, p=pico, etc.), as well as a few special
abbreviations like ``mV`` for millivolt, ``pF`` for picofarad, etc.

.. code:: ipython2

    1000*amp


::


    

    NameErrorTraceback (most recent call last)

    <ipython-input-4-1fbbf0aee6aa> in <module>()
    ----> 1 1000*amp
    

    NameError: name 'amp' is not defined


.. code:: ipython2

    1e6*volt


::


    

    NameErrorTraceback (most recent call last)

    <ipython-input-5-b9f47fe99cdc> in <module>()
    ----> 1 1e6*volt
    

    NameError: name 'volt' is not defined


.. code:: ipython2

    1000*namp


::


    

    NameErrorTraceback (most recent call last)

    <ipython-input-6-39a48d6121d2> in <module>()
    ----> 1 1000*namp
    

    NameError: name 'namp' is not defined


Also note that combinations of units with work as expected:

.. code:: ipython2

    10*nA*5*Mohm


::


    

    NameErrorTraceback (most recent call last)

    <ipython-input-7-7f4b07364399> in <module>()
    ----> 1 10*nA*5*Mohm
    

    NameError: name 'nA' is not defined


And if you try to do something wrong like adding amps and volts, what
happens?

.. code:: ipython2

    5*amp+10*volt


::


    

    NameErrorTraceback (most recent call last)

    <ipython-input-8-245c0c0332d1> in <module>()
    ----> 1 5*amp+10*volt
    

    NameError: name 'amp' is not defined


If you haven’t see an error message in Python before that can look a bit
overwhelming, but it’s actually quite simple and it’s important to know
how to read these because you’ll probably see them quite often.

You should start at the bottom and work up. The last line gives the
error type ``DimensionMismatchError`` along with a more specific message
(in this case, you were trying to add together two quantities with
different SI units, which is impossible).

Working upwards, each of the sections starts with a filename
(e.g. ``C:\Users\Dan\...``) with possibly the name of a function, and
then a few lines surrounding the line where the error occurred (which is
identified with an arrow).

The last of these sections shows the place in the function where the
error actually happened. The section above it shows the function that
called that function, and so on until the first section will be the
script that you actually run. This sequence of sections is called a
traceback, and is helpful in debugging.

If you see a traceback, what you want to do is start at the bottom and
scan up the sections until you find your own file because that’s most
likely where the problem is. (Of course, your code might be correct and
Brian may have a bug in which case, please let us know on the email
support list.)

A simple model
--------------

Let’s start by defining a simple neuron model. In Brian, all models are
defined by systems of differential equations. Here’s a simple example of
what that looks like:

.. code:: ipython2

    tau = 10*ms
    eqs = '''
    dv/dt = (1-v)/tau : 1
    '''


::


    

    NameErrorTraceback (most recent call last)

    <ipython-input-9-3914980972e6> in <module>()
    ----> 1 tau = 10*ms
          2 eqs = '''
          3 dv/dt = (1-v)/tau : 1
          4 '''


    NameError: name 'ms' is not defined


In Python, the notation ``'''`` is used to begin and end a multi-line
string. So the equations are just a string with one line per equation.
The equations are formatted with standard mathematical notation, with
one addition. At the end of a line you write ``: unit`` where ``unit``
is the SI unit of that variable. Note that this is not the unit of the
two sides of the equation (which would be ``1/second``), but the unit of
the *variable* defined by the equation, i.e. in this case :math:`v`.

Now let’s use this definition to create a neuron.

.. code:: ipython2

    G = NeuronGroup(1, eqs)


::


    

    NameErrorTraceback (most recent call last)

    <ipython-input-10-74b836f776d1> in <module>()
    ----> 1 G = NeuronGroup(1, eqs)
    

    NameError: name 'NeuronGroup' is not defined


In Brian, you only create groups of neurons, using the class
``NeuronGroup``. The first two arguments when you create one of these
objects are the number of neurons (in this case, 1) and the defining
differential equations.

Let’s see what happens if we didn’t put the variable ``tau`` in the
equation:

.. code:: ipython2

    eqs = '''
    dv/dt = 1-v : 1
    '''
    G = NeuronGroup(1, eqs)
    run(100*ms)


::


    

    NameErrorTraceback (most recent call last)

    <ipython-input-11-97ed109f5888> in <module>()
          2 dv/dt = 1-v : 1
          3 '''
    ----> 4 G = NeuronGroup(1, eqs)
          5 run(100*ms)


    NameError: name 'NeuronGroup' is not defined


An error is raised, but why? The reason is that the differential
equation is now dimensionally inconsistent. The left hand side ``dv/dt``
has units of ``1/second`` but the right hand side ``1-v`` is
dimensionless. People often find this behaviour of Brian confusing
because this sort of equation is very common in mathematics. However,
for quantities with physical dimensions it is incorrect because the
results would change depending on the unit you measured it in. For time,
if you measured it in seconds the same equation would behave differently
to how it would if you measured time in milliseconds. To avoid this, we
insist that you always specify dimensionally consistent equations.

Now let’s go back to the good equations and actually run the simulation.

.. code:: ipython2

    start_scope()
    
    tau = 10*ms
    eqs = '''
    dv/dt = (1-v)/tau : 1
    '''
    
    G = NeuronGroup(1, eqs)
    run(100*ms)


::


    

    NameErrorTraceback (most recent call last)

    <ipython-input-12-88b90c52ec07> in <module>()
    ----> 1 start_scope()
          2 
          3 tau = 10*ms
          4 eqs = '''
          5 dv/dt = (1-v)/tau : 1


    NameError: name 'start_scope' is not defined


First off, ignore that ``start_scope()`` at the top of the cell. You’ll
see that in each cell in this tutorial where we run a simulation. All it
does is make sure that any Brian objects created before the function is
called aren’t included in the next run of the simulation.

Secondly, you’ll see that there is an “INFO” message about not
specifying the numerical integration method. This is harmless and just
to let you know what method we chose, but we’ll fix it in the next cell
by specifying the method explicitly.

So, what has happened here? Well, the command ``run(100*ms)`` runs the
simulation for 100 ms. We can see that this has worked by printing the
value of the variable ``v`` before and after the simulation.

.. code:: ipython2

    start_scope()
    
    G = NeuronGroup(1, eqs, method='exact')
    print('Before v = %s' % G.v[0])
    run(100*ms)
    print('After v = %s' % G.v[0])


::


    

    NameErrorTraceback (most recent call last)

    <ipython-input-13-3e2d05959f69> in <module>()
    ----> 1 start_scope()
          2 
          3 G = NeuronGroup(1, eqs, method='exact')
          4 print('Before v = %s' % G.v[0])
          5 run(100*ms)


    NameError: name 'start_scope' is not defined


By default, all variables start with the value 0. Since the differential
equation is ``dv/dt=(1-v)/tau`` we would expect after a while that ``v``
would tend towards the value 1, which is just what we see. Specifically,
we’d expect ``v`` to have the value ``1-exp(-t/tau)``. Let’s see if
that’s right.

.. code:: ipython2

    print('Expected value of v = %s' % (1-exp(-100*ms/tau)))


::


    

    NameErrorTraceback (most recent call last)

    <ipython-input-14-8f1e09d2cba6> in <module>()
    ----> 1 print('Expected value of v = %s' % (1-exp(-100*ms/tau)))
    

    NameError: name 'exp' is not defined


Good news, the simulation gives the value we’d expect!

Now let’s take a look at a graph of how the variable ``v`` evolves over
time.

.. code:: ipython2

    start_scope()
    
    G = NeuronGroup(1, eqs, method='exact')
    M = StateMonitor(G, 'v', record=True)
    
    run(30*ms)
    
    plot(M.t/ms, M.v[0])
    xlabel('Time (ms)')
    ylabel('v');


::


    

    NameErrorTraceback (most recent call last)

    <ipython-input-15-a489b277bd78> in <module>()
    ----> 1 start_scope()
          2 
          3 G = NeuronGroup(1, eqs, method='exact')
          4 M = StateMonitor(G, 'v', record=True)
          5 


    NameError: name 'start_scope' is not defined


This time we only ran the simulation for 30 ms so that we can see the
behaviour better. It looks like it’s behaving as expected, but let’s
just check that analytically by plotting the expected behaviour on top.

.. code:: ipython2

    start_scope()
    
    G = NeuronGroup(1, eqs, method='exact')
    M = StateMonitor(G, 'v', record=0)
    
    run(30*ms)
    
    plot(M.t/ms, M.v[0], 'C0', label='Brian')
    plot(M.t/ms, 1-exp(-M.t/tau), 'C1--',label='Analytic')
    xlabel('Time (ms)')
    ylabel('v')
    legend();


::


    

    NameErrorTraceback (most recent call last)

    <ipython-input-16-b0ce8a03d476> in <module>()
    ----> 1 start_scope()
          2 
          3 G = NeuronGroup(1, eqs, method='exact')
          4 M = StateMonitor(G, 'v', record=0)
          5 


    NameError: name 'start_scope' is not defined


As you can see, the blue (Brian) and dashed orange (analytic solution)
lines coincide.

In this example, we used the object ``StateMonitor`` object. This is
used to record the values of a neuron variable while the simulation
runs. The first two arguments are the group to record from, and the
variable you want to record from. We also specify ``record=0``. This
means that we record all values for neuron 0. We have to specify which
neurons we want to record because in large simulations with many neurons
it usually uses up too much RAM to record the values of all neurons.

Now try modifying the equations and parameters and see what happens in
the cell below.

.. code:: ipython2

    start_scope()
    
    tau = 10*ms
    eqs = '''
    dv/dt = (sin(2*pi*100*Hz*t)-v)/tau : 1
    '''
    
    # Change to Euler method because exact integrator doesn't work here
    G = NeuronGroup(1, eqs, method='euler')
    M = StateMonitor(G, 'v', record=0)
    
    G.v = 5 # initial value
    
    run(60*ms)
    
    plot(M.t/ms, M.v[0])
    xlabel('Time (ms)')
    ylabel('v');


::


    

    NameErrorTraceback (most recent call last)

    <ipython-input-17-6915abd7dcbf> in <module>()
    ----> 1 start_scope()
          2 
          3 tau = 10*ms
          4 eqs = '''
          5 dv/dt = (sin(2*pi*100*Hz*t)-v)/tau : 1


    NameError: name 'start_scope' is not defined


Adding spikes
-------------

So far we haven’t done anything neuronal, just played around with
differential equations. Now let’s start adding spiking behaviour.

.. code:: ipython2

    start_scope()
    
    tau = 10*ms
    eqs = '''
    dv/dt = (1-v)/tau : 1
    '''
    
    G = NeuronGroup(1, eqs, threshold='v>0.8', reset='v = 0', method='exact')
    
    M = StateMonitor(G, 'v', record=0)
    run(50*ms)
    plot(M.t/ms, M.v[0])
    xlabel('Time (ms)')
    ylabel('v');


::


    

    NameErrorTraceback (most recent call last)

    <ipython-input-18-82651eb8971a> in <module>()
    ----> 1 start_scope()
          2 
          3 tau = 10*ms
          4 eqs = '''
          5 dv/dt = (1-v)/tau : 1


    NameError: name 'start_scope' is not defined


We’ve added two new keywords to the ``NeuronGroup`` declaration:
``threshold='v>0.8'`` and ``reset='v = 0'``. What this means is that
when ``v>0.8`` we fire a spike, and immediately reset ``v = 0`` after
the spike. We can put any expression and series of statements as these
strings.

As you can see, at the beginning the behaviour is the same as before
until ``v`` crosses the threshold ``v>0.8`` at which point you see it
reset to 0. You can’t see it in this figure, but internally Brian has
registered this event as a spike. Let’s have a look at that.

.. code:: ipython2

    start_scope()
    
    G = NeuronGroup(1, eqs, threshold='v>0.8', reset='v = 0', method='exact')
    
    spikemon = SpikeMonitor(G)
    
    run(50*ms)
    
    print('Spike times: %s' % spikemon.t[:])


::


    

    NameErrorTraceback (most recent call last)

    <ipython-input-19-64d6412551f0> in <module>()
    ----> 1 start_scope()
          2 
          3 G = NeuronGroup(1, eqs, threshold='v>0.8', reset='v = 0', method='exact')
          4 
          5 spikemon = SpikeMonitor(G)


    NameError: name 'start_scope' is not defined


The ``SpikeMonitor`` object takes the group whose spikes you want to
record as its argument and stores the spike times in the variable ``t``.
Let’s plot those spikes on top of the other figure to see that it’s
getting it right.

.. code:: ipython2

    start_scope()
    
    G = NeuronGroup(1, eqs, threshold='v>0.8', reset='v = 0', method='exact')
    
    statemon = StateMonitor(G, 'v', record=0)
    spikemon = SpikeMonitor(G)
    
    run(50*ms)
    
    plot(statemon.t/ms, statemon.v[0])
    for t in spikemon.t:
        axvline(t/ms, ls='--', c='C1', lw=3)
    xlabel('Time (ms)')
    ylabel('v');


::


    

    NameErrorTraceback (most recent call last)

    <ipython-input-20-1c92fa1f9d75> in <module>()
    ----> 1 start_scope()
          2 
          3 G = NeuronGroup(1, eqs, threshold='v>0.8', reset='v = 0', method='exact')
          4 
          5 statemon = StateMonitor(G, 'v', record=0)


    NameError: name 'start_scope' is not defined


Here we’ve used the ``axvline`` command from ``matplotlib`` to draw an
orange, dashed vertical line at the time of each spike recorded by the
``SpikeMonitor``.

Now try changing the strings for ``threshold`` and ``reset`` in the cell
above to see what happens.

Refractoriness
--------------

A common feature of neuron models is refractoriness. This means that
after the neuron fires a spike it becomes refractory for a certain
duration and cannot fire another spike until this period is over. Here’s
how we do that in Brian.

.. code:: ipython2

    start_scope()
    
    tau = 10*ms
    eqs = '''
    dv/dt = (1-v)/tau : 1 (unless refractory)
    '''
    
    G = NeuronGroup(1, eqs, threshold='v>0.8', reset='v = 0', refractory=5*ms, method='exact')
    
    statemon = StateMonitor(G, 'v', record=0)
    spikemon = SpikeMonitor(G)
    
    run(50*ms)
    
    plot(statemon.t/ms, statemon.v[0])
    for t in spikemon.t:
        axvline(t/ms, ls='--', c='C1', lw=3)
    xlabel('Time (ms)')
    ylabel('v');


::


    

    NameErrorTraceback (most recent call last)

    <ipython-input-21-d65e3bba0653> in <module>()
    ----> 1 start_scope()
          2 
          3 tau = 10*ms
          4 eqs = '''
          5 dv/dt = (1-v)/tau : 1 (unless refractory)


    NameError: name 'start_scope' is not defined


As you can see in this figure, after the first spike, ``v`` stays at 0
for around 5 ms before it resumes its normal behaviour. To do this,
we’ve done two things. Firstly, we’ve added the keyword
``refractory=5*ms`` to the ``NeuronGroup`` declaration. On its own, this
only means that the neuron cannot spike in this period (see below), but
doesn’t change how ``v`` behaves. In order to make ``v`` stay constant
during the refractory period, we have to add ``(unless refractory)`` to
the end of the definition of ``v`` in the differential equations. What
this means is that the differential equation determines the behaviour of
``v`` unless it’s refractory in which case it is switched off.

Here’s what would happen if we didn’t include ``(unless refractory)``.
Note that we’ve also decreased the value of ``tau`` and increased the
length of the refractory period to make the behaviour clearer.

.. code:: ipython2

    start_scope()
    
    tau = 5*ms
    eqs = '''
    dv/dt = (1-v)/tau : 1
    '''
    
    G = NeuronGroup(1, eqs, threshold='v>0.8', reset='v = 0', refractory=15*ms, method='exact')
    
    statemon = StateMonitor(G, 'v', record=0)
    spikemon = SpikeMonitor(G)
    
    run(50*ms)
    
    plot(statemon.t/ms, statemon.v[0])
    for t in spikemon.t:
        axvline(t/ms, ls='--', c='C1', lw=3)
    axhline(0.8, ls=':', c='C2', lw=3)
    xlabel('Time (ms)')
    ylabel('v')
    print("Spike times: %s" % spikemon.t[:])


::


    

    NameErrorTraceback (most recent call last)

    <ipython-input-22-b3bddd59f027> in <module>()
    ----> 1 start_scope()
          2 
          3 tau = 5*ms
          4 eqs = '''
          5 dv/dt = (1-v)/tau : 1


    NameError: name 'start_scope' is not defined


So what’s going on here? The behaviour for the first spike is the same:
``v`` rises to 0.8 and then the neuron fires a spike at time 8 ms before
immediately resetting to 0. Since the refractory period is now 15 ms
this means that the neuron won’t be able to spike again until time 8 +
15 = 23 ms. Immediately after the first spike, the value of ``v`` now
instantly starts to rise because we didn’t specify
``(unless refractory)`` in the definition of ``dv/dt``. However, once it
reaches the value 0.8 (the dashed green line) at time roughly 8 ms it
doesn’t fire a spike even though the threshold is ``v>0.8``. This is
because the neuron is still refractory until time 23 ms, at which point
it fires a spike.

Note that you can do more complicated and interesting things with
refractoriness. See the full documentation for more details about how it
works.

Multiple neurons
----------------

So far we’ve only been working with a single neuron. Let’s do something
interesting with multiple neurons.

.. code:: ipython2

    start_scope()
    
    N = 100
    tau = 10*ms
    eqs = '''
    dv/dt = (2-v)/tau : 1
    '''
    
    G = NeuronGroup(N, eqs, threshold='v>1', reset='v=0', method='exact')
    G.v = 'rand()'
    
    spikemon = SpikeMonitor(G)
    
    run(50*ms)
    
    plot(spikemon.t/ms, spikemon.i, '.k')
    xlabel('Time (ms)')
    ylabel('Neuron index');


::


    

    NameErrorTraceback (most recent call last)

    <ipython-input-23-71c26b49d858> in <module>()
    ----> 1 start_scope()
          2 
          3 N = 100
          4 tau = 10*ms
          5 eqs = '''


    NameError: name 'start_scope' is not defined


This shows a few changes. Firstly, we’ve got a new variable ``N``
determining the number of neurons. Secondly, we added the statement
``G.v = 'rand()'`` before the run. What this does is initialise each
neuron with a different uniform random value between 0 and 1. We’ve done
this just so each neuron will do something a bit different. The other
big change is how we plot the data in the end.

As well as the variable ``spikemon.t`` with the times of all the spikes,
we’ve also used the variable ``spikemon.i`` which gives the
corresponding neuron index for each spike, and plotted a single black
dot with time on the x-axis and neuron index on the y-value. This is the
standard “raster plot” used in neuroscience.

Parameters
----------

To make these multiple neurons do something more interesting, let’s
introduce per-neuron parameters that don’t have a differential equation
attached to them.

.. code:: ipython2

    start_scope()
    
    N = 100
    tau = 10*ms
    v0_max = 3.
    duration = 1000*ms
    
    eqs = '''
    dv/dt = (v0-v)/tau : 1 (unless refractory)
    v0 : 1
    '''
    
    G = NeuronGroup(N, eqs, threshold='v>1', reset='v=0', refractory=5*ms, method='exact')
    M = SpikeMonitor(G)
    
    G.v0 = 'i*v0_max/(N-1)'
    
    run(duration)
    
    figure(figsize=(12,4))
    subplot(121)
    plot(M.t/ms, M.i, '.k')
    xlabel('Time (ms)')
    ylabel('Neuron index')
    subplot(122)
    plot(G.v0, M.count/duration)
    xlabel('v0')
    ylabel('Firing rate (sp/s)');


::


    

    NameErrorTraceback (most recent call last)

    <ipython-input-24-f566f8d048d8> in <module>()
    ----> 1 start_scope()
          2 
          3 N = 100
          4 tau = 10*ms
          5 v0_max = 3.


    NameError: name 'start_scope' is not defined


The line ``v0 : 1`` declares a new per-neuron parameter ``v0`` with
units ``1`` (i.e. dimensionless).

The line ``G.v0 = 'i*v0_max/(N-1)'`` initialises the value of v0 for
each neuron varying from 0 up to ``v0_max``. The symbol ``i`` when it
appears in strings like this refers to the neuron index.

So in this example, we’re driving the neuron towards the value ``v0``
exponentially, but when ``v`` crosses ``v>1``, it fires a spike and
resets. The effect is that the rate at which it fires spikes will be
related to the value of ``v0``. For ``v0<1`` it will never fire a spike,
and as ``v0`` gets larger it will fire spikes at a higher rate. The
right hand plot shows the firing rate as a function of the value of
``v0``. This is the I-f curve of this neuron model.

Note that in the plot we’ve used the ``count`` variable of the
``SpikeMonitor``: this is an array of the number of spikes each neuron
in the group fired. Dividing this by the duration of the run gives the
firing rate.

Stochastic neurons
------------------

Often when making models of neurons, we include a random element to
model the effect of various forms of neural noise. In Brian, we can do
this by using the symbol ``xi`` in differential equations. Strictly
speaking, this symbol is a “stochastic differential” but you can sort of
thinking of it as just a Gaussian random variable with mean 0 and
standard deviation 1. We do have to take into account the way stochastic
differentials scale with time, which is why we multiply it by
``tau**-0.5`` in the equations below (see a textbook on stochastic
differential equations for more details). Note that we also changed the
``method`` keyword argument to use ``'euler'`` (which stands for the
`Euler-Maruyama
method <https://en.wikipedia.org/wiki/Euler%E2%80%93Maruyama_method>`__);
the ``'exact'`` method that we used earlier is not applicable to
stochastic differential equations.

.. code:: ipython2

    start_scope()
    
    N = 100
    tau = 10*ms
    v0_max = 3.
    duration = 1000*ms
    sigma = 0.2
    
    eqs = '''
    dv/dt = (v0-v)/tau+sigma*xi*tau**-0.5 : 1 (unless refractory)
    v0 : 1
    '''
    
    G = NeuronGroup(N, eqs, threshold='v>1', reset='v=0', refractory=5*ms, method='euler')
    M = SpikeMonitor(G)
    
    G.v0 = 'i*v0_max/(N-1)'
    
    run(duration)
    
    figure(figsize=(12,4))
    subplot(121)
    plot(M.t/ms, M.i, '.k')
    xlabel('Time (ms)')
    ylabel('Neuron index')
    subplot(122)
    plot(G.v0, M.count/duration)
    xlabel('v0')
    ylabel('Firing rate (sp/s)');


::


    

    NameErrorTraceback (most recent call last)

    <ipython-input-25-e38c73ce4585> in <module>()
    ----> 1 start_scope()
          2 
          3 N = 100
          4 tau = 10*ms
          5 v0_max = 3.


    NameError: name 'start_scope' is not defined


That’s the same figure as in the previous section but with some noise
added. Note how the curve has changed shape: instead of a sharp jump
from firing at rate 0 to firing at a positive rate, it now increases in
a sigmoidal fashion. This is because no matter how small the driving
force the randomness may cause it to fire a spike.

End of tutorial
---------------

That’s the end of this part of the tutorial. The cell below has another
example. See if you can work out what it is doing and why. Try adding a
``StateMonitor`` to record the values of the variables for one of the
neurons to help you understand it.

You could also try out the things you’ve learned in this cell.

Once you’re done with that you can move on to the next tutorial on
Synapses.

.. code:: ipython2

    start_scope()
    
    N = 1000
    tau = 10*ms
    vr = -70*mV
    vt0 = -50*mV
    delta_vt0 = 5*mV
    tau_t = 100*ms
    sigma = 0.5*(vt0-vr)
    v_drive = 2*(vt0-vr)
    duration = 100*ms
    
    eqs = '''
    dv/dt = (v_drive+vr-v)/tau + sigma*xi*tau**-0.5 : volt
    dvt/dt = (vt0-vt)/tau_t : volt
    '''
    
    reset = '''
    v = vr
    vt += delta_vt0
    '''
    
    G = NeuronGroup(N, eqs, threshold='v>vt', reset=reset, refractory=5*ms, method='euler')
    spikemon = SpikeMonitor(G)
    
    G.v = 'rand()*(vt0-vr)+vr'
    G.vt = vt0
    
    run(duration)
    
    _ = hist(spikemon.t/ms, 100, histtype='stepfilled', facecolor='k', weights=list(ones(len(spikemon))/(N*defaultclock.dt)))
    xlabel('Time (ms)')
    ylabel('Instantaneous firing rate (sp/s)');


::


    

    NameErrorTraceback (most recent call last)

    <ipython-input-26-ed364f7aa4b4> in <module>()
    ----> 1 start_scope()
          2 
          3 N = 1000
          4 tau = 10*ms
          5 vr = -70*mV


    NameError: name 'start_scope' is not defined

