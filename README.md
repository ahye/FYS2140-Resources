FYS2140-Resources
=================

Here you can find all python resources for the FYS2140 quantum physics
course. All scripts are kept in the `src`-folder.

## How to download
In order to download the files you can on any linux-machine type the
following command in the terminal:

```bash
$ git clone https://github.com/ahye/FYS2140-Resources.git
```

This will create a new folder `FYS2140-Resources` where everything is
stored. You can also check back to this website in order to watch for
changes or updates to the file. The right column in the file tree says
when the file was last updated. A short message describing the changes
is also included.

You can also use the button to the right labeled "Download ZIP". This
leaves you with a ZIP-archive that you'll need to extract.

## Necessary software
The scripts using the Crank Nicolson algorithm are dependent upon the
`scipy` library.

### On UiO computers
You can use the computers in the computer lab in room FV329 in the
physics building. You can also use `ssh` to login on those computers,
but animation will go slowly.

```bash
$ ssh -XY <username>@spole.uio.no
```

Where `spole` is the name of one of the computers. This can be
replaced by any of the others. You can find other names by looking at
the white labels on the computers.

### On your own computer
Usually `scipy` is already installed if you use a mac or a linux based
computer. If not, you can on linux machines write the following
command in the terminal:

```bash
$ sudo apt-get install python-scipy
```

On a mac you can install the newest Enthought Python distribution.  This is free
(32-bit) and can be found [here](https://www.enthought.com/products/epd/free/).
If you want the 64-bit version you can ask your group teacher --- the university
has licenses you can use.

## Contribute!
The goal of this package is to show some examples on how the
Schrodinger equation can be solved using different numerical methods.
What the different methods do better than other and so on.

If you have any ideas on how this can be made better. Please do
contact me by e-mail at `benedicte[at]ifi[dot]uio[dot]no` or send
pull requests with describing messages.
