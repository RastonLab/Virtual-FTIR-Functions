from app import app
import sys

if __name__ == "__main__":

    debug = False
    if len(sys.argv) > 1 and sys.argv[1].lower() == "debug":
        debug = bool(sys.argv[1])
    
    app.run(debug=debug)
