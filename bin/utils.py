import os
import re
import shutil
import subprocess

def get_files_path(folder_path):
    """Get the absolute path of all files in the current directory"""
    file_paths = []
    for item in os.listdir(folder_path):
        item_path = os.path.join(folder_path, item)
        if os.path.isfile(item_path):
            absolute_path = os.path.abspath(item_path)
            file_paths.append(absolute_path)
    return file_paths


def get_files_prefix(file_paths, file_extension=None):
    """
    Get the prefix of all files in the input list
    
    file_extension: specifies the file extension (such as '.bam', '.txt', etc.). If None, all files will be processed.
    """
    file_prefixes = []
    none_paths = []
    for path in file_paths:
        if os.path.isfile(path):
            filename = os.path.basename(path)
            if file_extension is not None and not filename.endswith(file_extension):
                file_prefixes.append(None)
                none_paths.append(path)
            else:
                prefix = os.path.splitext(filename)[0]
                file_prefixes.append(prefix)

    if None in file_prefixes:
        print(f'Error: Some files do not match the specified file extension "{file_extension}" and will be skipped.')
        print(f'Files that do not match: {none_paths}')
        return None

    if len(set(file_prefixes))!=len(file_prefixes):
        duplicate_info = []
        for i, prefix in enumerate(file_prefixes):
            duplicate_indices = [j for j, p in enumerate(file_prefixes) if p == prefix]
            if len(duplicate_indices) > 1:
                duplicate_paths = [file_paths[j] for j in duplicate_indices]
                info = f"'{prefix}': {duplicate_paths}"
                if info not in duplicate_info:
                    duplicate_info.append(info)
        print('Error: Duplicate file prefixes detected.')
        print(f'Duplicate prefixes and their corresponding paths: {", ".join(duplicate_info)}')
        return None
    
    return file_prefixes


def get_merged_combs(unique_combs) -> list:
    merged_combs_raw = []
    for unique_comb_i in unique_combs:
        target1, target2 = unique_comb_i.split("-")
        if target1 <= target2:
            merged_combs_raw.append(f"{target1}-{target2}")
    merged_combs = list(set(merged_combs_raw))  # Remove duplicates
    print(f"Unique combinations: {merged_combs}")
    return merged_combs


def get_fastq_prefix(file_paths, ext=".fastq.gz"):
    PAIRED_END_PATTERNS = [
        r'_R[12]$',    # _R1, _R2
        r'\.R[12]$',   # .R1, .R2
        r'_L00[12]_R[12]$',      # _L001_R1, _L001_R2
        r'_L00[12]_[12]$',       # _L001_1, _L001_2
        r'\.L00[12]\.[12]$',     # .L001.1, .L001.2
    ]
    FORWARD_PATTERNS = [
        r'_R1$',    
        r'\.R1$',   
        r'_L001_R1$',      
        r'_L001_1$',       
        r'\.L001\.1$',     
    ]
    
    if len(file_paths) % 2 != 0:
        print(f"ERROR: Expected even number of files for pairing, got {len(file_paths)}")
        return None
    expected_pairs = len(file_paths) // 2

    prefix_groups = {}
    unmatched_files = []
    pattern_matches = {}
    for path in file_paths:
        basename = os.path.basename(path).replace(ext, '') 
        cleaned = basename
        matched = False
        matched_pattern_index = None
        # Try each pattern and remove the first match
        for idx, pattern in enumerate(PAIRED_END_PATTERNS):
            if re.search(pattern, cleaned):
                cleaned = re.sub(pattern, '', cleaned)
                if pattern not in pattern_matches:
                    pattern_matches[pattern] = []  # Track pattern usage for debugging
                pattern_matches[pattern].append(basename)
                matched = True
                matched_pattern_index = idx
                break
        if not matched:
            unmatched_files.append(basename)
        if cleaned not in prefix_groups:
            prefix_groups[cleaned] = []
        prefix_groups[cleaned].append((path, basename, matched_pattern_index))

    if unmatched_files:
        print(f"ERROR: {len(unmatched_files)} files did not match any pattern:")
        for unmatched in unmatched_files:
             print(f"  - {unmatched}")
        return None
    
    if len(pattern_matches) == 1:
        pattern_name = list(pattern_matches.keys())[0]
        print(f"All files matched pattern: {pattern_name}")
    elif len(pattern_matches) > 1:
        print("WARNING: Multiple patterns detected!")
        for pattern, matches in pattern_matches.items():
            print(f"  Pattern '{pattern}' matched {len(matches)} files:")
            for match in matches:
                print(f"    - {match}")
    
    if len(prefix_groups) != expected_pairs:
        print(f"ERROR: Expected {expected_pairs} unique prefixes, got {len(prefix_groups)}")
        print("Prefix groups found:")
        for prefix, files in prefix_groups.items():
            print(f"  {prefix}: {len(files)} files")
        return None
    
    result = {}
    for prefix, file_info_list in prefix_groups.items():
        if len(file_info_list) != 2:
            print(f"ERROR: Prefix '{prefix}' has {len(file_info_list)} files, expected exactly 2")
            return None
        fwd_file, rev_file = None, None
        for path, basename, pattern_index in file_info_list:
            if re.search(FORWARD_PATTERNS[pattern_index], basename):
                if fwd_file is not None:
                    print(f"ERROR: Multiple forward reads found for prefix '{prefix}'")
                    return None
                fwd_file = path
            else:
                if rev_file is not None:
                    print(f"ERROR: Multiple reverse reads found for prefix '{prefix}'")
                    return None
                rev_file = path
        if fwd_file is None or rev_file is None:
            print(f"ERROR: Could not identify forward/reverse pair for prefix '{prefix}'")
            return None
        result[prefix] = {'fwd': fwd_file, 'rev': rev_file}
    
    print(f"Successfully paired {len(result)} sample pairs")
    return result


def get_picard_jar_path():
    """Capture the jar path from current conda environment"""
    picard_script = shutil.which('picard')
    if not picard_script or not os.path.exists(picard_script): 
        print("Error: picard is not installed or not in PATH")
        return None
    
    source = picard_script
    while os.path.islink(source):
        current_dir = os.path.dirname(os.path.abspath(source))
        link_target = os.readlink(source)
        if not os.path.isabs(link_target):
            source = os.path.join(current_dir, link_target)
        else:
            source = link_target
    jar_dir = os.path.dirname(os.path.abspath(source))
    jar_path = os.path.join(jar_dir, 'picard.jar')
    if os.path.exists(jar_path):
        return jar_path
    else:
        print(f"Error: jar file does not exist in the expected location {jar_path}")
        return None
