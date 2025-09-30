# Python script to diff two UM job outputs and check the namelist inputs
import argparse
import re
import difflib
import sys

def filter_file_content(file_path: str, pattern: str) -> list[str]:
    """
    Reads a file and filters its lines based on a regex pattern (simulating grep).

    Args:
        file_path: The path to the input text file.
        pattern: The regular expression pattern to match against each line.

    Returns:
        A list of lines from the file that match the pattern, stripped of leading/trailing whitespace.
        Returns an empty list if the file cannot be read.
    """
    matching_lines = []
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            print(f"-> Reading and filtering content from {file_path}...")
            # Compile the regex pattern for efficiency
            compiled_pattern = re.compile(pattern)
            
            for line in f:
                # Use search to find the pattern anywhere in the line, 
                # similar to standard grep behavior
                if compiled_pattern.search(line):
                    # We strip whitespace for cleaner diff output, but keep the newline
                    # only if you want the diff to include line breaks in the comparison.
                    # For a line-by-line comparison where lines might not have trailing
                    # spaces, it's often best to keep the newline character 
                    # for accurate difflib output, or remove it entirely.
                    # Here, we remove it to compare the core content only:
                    matching_lines.append(line.rstrip())

    except FileNotFoundError:
        print(f"Error: File not found at {file_path}", file=sys.stderr)
    except Exception as e:
        print(f"An unexpected error occurred while processing {file_path}: {e}", file=sys.stderr)
    
    return matching_lines


def generate_subset_diff(file1_path: str, file2_path: str, pattern: str):
    """
    Filters two files based on a pattern and prints the unified difference report.
    """
    
    # 1. Filter content from both files
    content1 = filter_file_content(file1_path, pattern)
    content2 = filter_file_content(file2_path, pattern)

    if not content1 and not content2:
        print("\nNo matching lines found in either file, or file errors occurred. Exiting.")
        return
        
    print(f"\n--- Subset 1 ({file1_path}) has {len(content1)} matching lines.")
    print(f"--- Subset 2 ({file2_path}) has {len(content2)} matching lines.")
    print(f"--- Generating unified diff using pattern: '{pattern}'\n")

    # 2. Calculate the difference using difflib
    # The 'fromfile' and 'tofile' arguments provide headers for the diff output
    diff = difflib.unified_diff(
        content1, 
        content2, 
        fromfile=f'a/{file1_path} (filtered by: {pattern})', 
        tofile=f'b/{file2_path} (filtered by: {pattern})',
        lineterm='' # Prevent adding extra newlines since we removed them during filtering
    )

    # 3. Print the diff report
    for line in diff:
        print(line)

    # 4. Is the content the same?
    if sorted(content2) == sorted(content1):
        print (f" Lines flagged by '{pattern}' are identical")

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="Generate a unified diff between subsets of two files, filtered by a regular expression pattern (simulating grep)."
    )
    parser.add_argument('file1', help="Path to the first text file (File A).")
    parser.add_argument('file2', help="Path to the second text file (File B).")
    parser.add_argument('pattern', help="The regular expression pattern to filter lines in both files.")
    
    args = parser.parse_args()
    
    # Run the difference generation
    generate_subset_diff(args.file1, args.file2, args.pattern)

    # Sort and find missing and matching patterns without whitespace

# Example Usage (after saving as subset_diff.py):
# python subset_diff.py file_a.txt file_b.txt 'ERROR|WARNING'

