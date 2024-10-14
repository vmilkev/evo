import evoped

# ----------------------------------
def example_Amat():

    amat = evoped.Amat()
    help(amat)

# ----------------------------------
def example_Gmat():

    gmat = evoped.Gmat()
    help(gmat)

# ----------------------------------
def example_Hmat():

    hmat = evoped.Hmat()
    help(hmat)

# ----------------------------------    
def main():

    example_Amat()
    print("completed: example_Amat")

    example_Gmat()
    print("completed: example_Gmat")

    example_Hmat()
    print("completed: example_Hmat")

# ----------------------------------	
if __name__ == '__main__':
    main()