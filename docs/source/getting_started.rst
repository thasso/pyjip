Getting startet guide
=====================

Hello world
-----------
In order to create your first jip script, create a file `hello_world.jip` with
the following content::
    
    #!/usr/bin/env jip
    echo "hello world"

By default, jip script commands are interpreted by bash. Make the file executable
and you can run it::

    $> chmod +x hello_world.jip
    $> ./hello_world.jip
    hello world
    $>

