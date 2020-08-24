The R console itself can be used as a fancy calculator. For example:

    1 + 2

    ## [1] 3

    2 * 6

    ## [1] 12

    5 / 5

    ## [1] 1

    2^2

    ## [1] 4

    log2(2)

    ## [1] 1

An important unit in programming is the *variable*. Variables hold
different *objects*, such as numbers, letters, and even more complex
types of data. To assign a variable in R, you follow the format
`variable_name <- object`.

Here I am storing a number as a variable.

    favorite_number <- 12

For all intents and purposes, `favorite_number` is now the same as the
number 12.

    favorite_number + 3

    ## [1] 15

As I eluded to earlier, there are different object types in R. Some of
the basic object types include integer, numeric, character, and factor
(which will be talked about later). You can check what type of object
something is by using the `class` function.

    # Integers are whole numbers, and are the number followed by the letter L.
    class(12L)

    ## [1] "integer"

    # Numeric is any number without an L.
    class(12)

    ## [1] "numeric"

    class(1.5)

    ## [1] "numeric"

    # Characters are letters.
    class("Bob")

    ## [1] "character"

    # Characters can also be numbers if you wrap them in quotes.
    class("1")

    ## [1] "character"
