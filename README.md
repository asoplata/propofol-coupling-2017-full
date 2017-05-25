# (Soplata et al., 2017) TC, RE, and synaptic DynaSim full simulation reproduction files.

This contains the complete [DynaSim](https://github.com/DynaSim/DynaSim) code
files (in `dynasim`), model mechanism files (in
`models/personal/propofol-coupling-2017-mechanisms`), and simulation runscripts
(in `runscripts`) needed for simulation of the thalamus of (Soplata et al.,
2017), inherited from (Ching et al., 2010). Note that you should NOT run these
scripts unless you have access to a parallel computing cluster that DynaSim
integrates with (and all the MATLAB licenses that that entails)! The main batch
of these involve almost 10,000 simulations of ~100 relatively complex cells (and
their combinatorial synaptic connections) for almost 10 seconds at a high
fixed-step resolution (Euler, 0.01 ms), and each individual simulation creates
over 1GB of data! All simulations were set to run with a default memory/RAM
allowance of 256GB. You will NOT be able to run all these simulations on your
personal computer in a reasonable amount of time, unless you are a time traveler
from the year 2054.

Adding the `dynasim` directory and all its subdirectories to your MATLAB path
should enable you to run all the simulations for the
computational thalamus from:

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
other words, the code contained here, used in (Soplata et al., 2017), is
**more** correct than the equations in the Supplementary Information of (Ching
et al., 2010).

Also note that, while this will create all of the subplots used in the figures
in the paper INDIVIDUALLY, the subplots were reformatted in Inkscape by hand. No
data or coloration was changed, and axis tick values were replaced with the same
values, but formatted differently to look better and for a more consistent
style.

## Original reproduction

These mechanism files are approximately as they were used originally, which
means this copy of DynaSim is using commit
[4a20467](https://github.com/DynaSim/DynaSim/commit/4a20467848a82673492ee06322acd3505e8c1788);
however, the version used for each batch of simulation may have varied, and so
feel free to contact me if there are any errors preventing simulations from
running from scratch.

If you want just the mechanism files alone, go to
[propofol-coupling-2017-mechanisms](https://github.com/asoplata/propofol-coupling-2017-mechanisms).
