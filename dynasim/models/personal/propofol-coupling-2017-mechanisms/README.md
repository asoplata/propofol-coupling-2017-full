# (Soplata et al., 2017) TC, RE, and synaptic DynaSim mechanism files.

 DynaSim-compatible mechanism files for simulation of the thalamus of (Soplata
 et al., 2017), inherited from (Ching et al., 2010).

Adding these mechanism files and associated functions into where you keep your mechanism files for [DynaSim](https://github.com/DynaSim/DynaSim), e.g. `/your/path/to/dynasim/models`, should enable you to simulate the computational thalamus from:

    Soplata, A. E., McCarthy, M. M., Sherfey, J. Lee, S., Purdon, P. L.,
    Brown, E. N., & Kopell, N. J. (2017). Thalamocortical control of propofol
    phase-amplitude coupling. In preparation.

This model was inherited from:

    Ching, S., Cimenser, A., Purdon, P. L., Brown, E. N., & Kopell, N. J.
    (2010). Thalamocortical model for a propofol-induced alpha-rhythm
    associated with loss of consciousness. Proceedings of the National
    Academy of Sciences, 107(52), 22665–22670.
    http://doi.org/10.1073/pnas.1017069108

Note that this code diverges from the given equations of (Ching et al., 2010)
due to typos and errors in the original equations. The code contained here has
been ground-truthed and cross-checked across both the original code run for the
paper and the source material that THAT code is based on, several times. In
other words, the code contained here is **more** correct than the equations in
the paper.

## Install

The easiest way to get started with this is

1. Install DynaSim ([see here for
   instructions](https://github.com/DynaSim/DynaSim/wiki/Installation))
2. `git clone` this repo into '/your/path/to/dynasim/models', i.e. the 'models'
   subdirectory of your copy of the DynaSim repo.
3. Set your own data directory inside the sample runscript `tcre_run_script.m`,
   make your own changes, etc.
4. Believe it or not...that should be it! You should be able to start MATLAB in
   any directory and run this script successfully! Let me know if there are
   problems, at austin.soplata 'at symbol' gmail 'dot' com

Note that there are extra, unused synaptic mechanism files that do things like
simulate a persistent AMPAergic spike train down to the cells, etc., in order to
mimic the cortex's input.

## Original reproduction

These mechanism files are AS THEY WERE USED originally, which means they were
used under a version of DynaSim likely using the following commit:
[4a20467](https://github.com/DynaSim/DynaSim/commit/4a20467848a82673492ee06322acd3505e8c1788)

If you want a COMPLETE installation to reproduce the simulations used in the
paper, go to
[propofol-coupling-2017-full](https://github.com/asoplata/propofol-coupling-2017-full).
