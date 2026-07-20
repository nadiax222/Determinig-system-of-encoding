from pathlib import Path


def read_data_from_file(file_path):
    """
    Reads a DNA sequence and a quality string from a text file.

    Expected format:
    line 1 - DNA sequence
    line 2 - encoded quality string
    """
    file_path = Path(file_path)

    if not file_path.exists():
        raise FileNotFoundError(
            f"The file '{file_path}' does not exist."
        )

    if not file_path.is_file():
        raise ValueError(
            f"The path '{file_path}' is not a file."
        )

    try:
        with open(file_path, "r", encoding="utf-8") as file:
            lines = file.read().splitlines()

    except OSError as error:
        raise OSError(
            f"An error occurred while reading the file: {error}"
        )

    non_empty_lines = [
        line.strip()
        for line in lines
        if line.strip()
    ]

    if len(non_empty_lines) < 2:
        raise ValueError(
            "The file must contain at least two non-empty lines:\n"
            "line 1 - DNA sequence\n"
            "line 2 - quality string"
        )

    sequence = non_empty_lines[0].replace(" ", "").upper()
    quality_string = non_empty_lines[1]

    return sequence, quality_string


def validate_data(sequence, quality_string):
    """
    Validates the DNA sequence and quality string.

    Different lengths are allowed, but only matching positions
    will be analysed.
    """
    if not sequence:
        raise ValueError(
            "The DNA sequence is empty."
        )

    if not quality_string:
        raise ValueError(
            "The quality string is empty."
        )

    valid_bases = set("ACGTN")

    invalid_bases = [
        base
        for base in sequence
        if base not in valid_bases
    ]

    if invalid_bases:
        raise ValueError(
            "The DNA sequence contains invalid bases: "
            + ", ".join(sorted(set(invalid_bases)))
        )

    if len(sequence) != len(quality_string):
        print(
            "\nWarning: The DNA sequence and quality string "
            "have different lengths."
        )
        print(f"Sequence length: {len(sequence)}")
        print(f"Quality string length: {len(quality_string)}")
        print(
            "Only positions containing both a nucleotide "
            "and a quality character will be analysed."
        )


def detect_encoding(quality_string):
    """
    Attempts to detect Phred+33 or Phred+64 encoding.

    The ASCII ranges overlap, so some results may be ambiguous.
    """
    ascii_values = [
        ord(character)
        for character in quality_string
    ]

    min_ascii = min(ascii_values)
    max_ascii = max(ascii_values)

    if min_ascii < 59:
        encoding = "Phred+33"

    elif min_ascii >= 64:
        encoding = "Phred+64 or ambiguous"

    else:
        encoding = "Ambiguous, probably Phred+33"

    return encoding, min_ascii, max_ascii


def quality_to_phred(quality_string, encoding):
    """
    Converts quality characters into numerical Phred scores.
    """
    if encoding == "Phred+33":
        offset = 33

    elif encoding == "Phred+64":
        offset = 64

    else:
        raise ValueError(
            "Encoding must be 'Phred+33' or 'Phred+64'."
        )

    return [
        ord(character) - offset
        for character in quality_string
    ]


def phred_to_quality_string(phred_scores, encoding):
    """
    Converts numerical Phred scores into an encoded quality string.
    """
    if encoding == "Phred+33":
        offset = 33

    elif encoding == "Phred+64":
        offset = 64

    else:
        raise ValueError(
            "Encoding must be 'Phred+33' or 'Phred+64'."
        )

    return "".join(
        chr(score + offset)
        for score in phred_scores
    )


def calculate_quality_statistics(phred_scores):
    """
    Calculates basic Phred quality statistics.
    """
    if not phred_scores:
        raise ValueError(
            "The list of Phred scores is empty."
        )

    minimum_score = min(phred_scores)
    maximum_score = max(phred_scores)
    average_score = sum(phred_scores) / len(phred_scores)

    q20_count = sum(
        score >= 20
        for score in phred_scores
    )

    q30_count = sum(
        score >= 30
        for score in phred_scores
    )

    q20_percentage = (
        q20_count / len(phred_scores)
    ) * 100

    q30_percentage = (
        q30_count / len(phred_scores)
    ) * 100

    return {
        "minimum": minimum_score,
        "maximum": maximum_score,
        "average": average_score,
        "q20_percent": q20_percentage,
        "q30_percent": q30_percentage
    }


def find_low_quality_positions(phred_scores, threshold=20):
    """
    Returns positions with Phred score below the threshold.

    Positions are numbered starting from 1.
    """
    return [
        position
        for position, score in enumerate(
            phred_scores,
            start=1
        )
        if score < threshold
    ]


def find_low_quality_regions(phred_scores, threshold=20):
    """
    Finds continuous low-quality regions.

    Positions are numbered starting from 1.
    """
    low_quality_regions = []
    start = None

    for position, score in enumerate(
        phred_scores,
        start=1
    ):
        if score < threshold:
            if start is None:
                start = position

        elif start is not None:
            low_quality_regions.append(
                (start, position - 1)
            )
            start = None

    if start is not None:
        low_quality_regions.append(
            (start, len(phred_scores))
        )

    return low_quality_regions


def create_quality_table(
    sequence,
    quality_string,
    phred_scores
):
    """
    Creates a table containing sequence and quality information.
    """
    table_lines = [
        "Position\tBase\tQuality character\tPhred score"
    ]

    for position, (
        base,
        quality_character,
        phred_score
    ) in enumerate(
        zip(
            sequence,
            quality_string,
            phred_scores
        ),
        start=1
    ):
        table_lines.append(
            f"{position}\t"
            f"{base}\t"
            f"{quality_character}\t"
            f"{phred_score}"
        )

    return "\n".join(table_lines)


def trim_low_quality_ends(
    sequence,
    quality_string,
    phred_scores,
    threshold=20
):
    """
    Removes low-quality bases from both ends.

    Internal low-quality positions are not removed.
    """
    start = 0
    end = len(phred_scores)

    while (
        start < end
        and phred_scores[start] < threshold
    ):
        start += 1

    while (
        end > start
        and phred_scores[end - 1] < threshold
    ):
        end -= 1

    trimmed_sequence = sequence[start:end]
    trimmed_quality_string = quality_string[start:end]
    trimmed_phred_scores = phred_scores[start:end]

    return (
        trimmed_sequence,
        trimmed_quality_string,
        trimmed_phred_scores
    )


def create_report(
    filename,
    original_sequence_length,
    original_quality_length,
    sequence,
    quality_string,
    selected_encoding,
    min_ascii,
    max_ascii,
    statistics,
    low_quality_positions,
    low_quality_regions,
    trimmed_sequence,
    trimmed_quality_string,
    trimmed_scores,
    phred_scores
):
    """
    Creates the complete text report.
    """
    report_lines = [
        "=" * 65,
        "DNA QUALITY ANALYSIS REPORT",
        "=" * 65,
        f"Input file: {filename}",
        f"Original DNA sequence length: {original_sequence_length}",
        f"Original quality string length: {original_quality_length}",
        f"Analysed positions: {len(sequence)}",
        f"Analysed sequence: {sequence}",
        f"Analysed quality string: {quality_string}",
        f"Selected encoding: {selected_encoding}",
        f"Minimum ASCII value: {min_ascii}",
        f"Maximum ASCII value: {max_ascii}",
        "",
        "QUALITY STATISTICS",
        "-" * 65,
        f"Average Phred score: {statistics['average']:.2f}",
        f"Minimum Phred score: {statistics['minimum']}",
        f"Maximum Phred score: {statistics['maximum']}",
        (
            "Bases with quality Q20 or higher: "
            f"{statistics['q20_percent']:.2f}%"
        ),
        (
            "Bases with quality Q30 or higher: "
            f"{statistics['q30_percent']:.2f}%"
        ),
        "",
        "LOW-QUALITY POSITIONS",
        "-" * 65
    ]

    if low_quality_positions:
        report_lines.append(
            ", ".join(
                str(position)
                for position in low_quality_positions
            )
        )
    else:
        report_lines.append(
            "No positions below Q20 were found."
        )

    report_lines.extend([
        "",
        "LOW-QUALITY REGIONS",
        "-" * 65
    ])

    if low_quality_regions:
        for start, end in low_quality_regions:
            if start == end:
                report_lines.append(
                    f"Position {start}"
                )
            else:
                report_lines.append(
                    f"Positions {start}-{end}"
                )
    else:
        report_lines.append(
            "No regions below Q20 were found."
        )

    report_lines.extend([
        "",
        "TRIMMING RESULTS",
        "-" * 65,
        f"Analysed sequence length: {len(sequence)}",
        f"Trimmed sequence length: {len(trimmed_sequence)}",
        (
            f"Trimmed sequence: "
            f"{trimmed_sequence or 'Empty'}"
        ),
        (
            f"Trimmed quality string: "
            f"{trimmed_quality_string or 'Empty'}"
        ),
        (
            "Trimmed Phred scores: "
            + (
                ", ".join(
                    str(score)
                    for score in trimmed_scores
                )
                if trimmed_scores
                else "None"
            )
        ),
        "",
        "DETAILED QUALITY TABLE",
        "-" * 65,
        create_quality_table(
            sequence,
            quality_string,
            phred_scores
        ),
        "=" * 65
    ])

    return "\n".join(report_lines)


def save_report(report, filename):
    """
    Saves the analysis report to a text file.
    """
    if not filename:
        filename = "quality_report.txt"

    if not filename.lower().endswith(".txt"):
        filename += ".txt"

    output_directory = Path("reports")
    output_directory.mkdir(exist_ok=True)

    output_path = output_directory / filename

    try:
        with open(
            output_path,
            "w",
            encoding="utf-8"
        ) as file:
            file.write(report)

        print(
            f"\nReport saved to:\n{output_path.resolve()}"
        )

    except OSError as error:
        print(
            f"\nThe report could not be saved: {error}"
        )


def main():
    try:
        filename = input(
            "Enter the filename or full path: "
        ).strip()

        sequence, quality_string = read_data_from_file(
            filename
        )

        validate_data(
            sequence,
            quality_string
        )

        original_sequence_length = len(sequence)
        original_quality_length = len(quality_string)

        analysis_length = min(
            original_sequence_length,
            original_quality_length
        )

        sequence = sequence[:analysis_length]
        quality_string = quality_string[:analysis_length]

        print("\nData read successfully.")
        print(
            f"Original sequence length: "
            f"{original_sequence_length}"
        )
        print(
            f"Original quality string length: "
            f"{original_quality_length}"
        )
        print(
            f"Positions used in analysis: "
            f"{analysis_length}"
        )
        print(f"Analysed sequence: {sequence}")
        print(f"Analysed quality string: {quality_string}")

        detected_encoding, min_ascii, max_ascii = (
            detect_encoding(quality_string)
        )

        print(
            f"\nDetected encoding: "
            f"{detected_encoding}"
        )
        print(f"Minimum ASCII value: {min_ascii}")
        print(f"Maximum ASCII value: {max_ascii}")

        if detected_encoding == "Phred+33":
            selected_encoding = "Phred+33"

        elif detected_encoding == "Phred+64 or ambiguous":
            choice = input(
                "\nThe encoding may be ambiguous.\n"
                "Enter 33 for Phred+33 or 64 for Phred+64: "
            ).strip()

            if choice == "64":
                selected_encoding = "Phred+64"
            else:
                selected_encoding = "Phred+33"

        else:
            selected_encoding = "Phred+33"

        phred_scores = quality_to_phred(
            quality_string,
            selected_encoding
        )

        if any(score < 0 for score in phred_scores):
            raise ValueError(
                "Negative Phred scores were detected. "
                "The selected encoding is probably incorrect."
            )

        statistics = calculate_quality_statistics(
            phred_scores
        )

        low_quality_positions = find_low_quality_positions(
            phred_scores,
            threshold=20
        )

        low_quality_regions = find_low_quality_regions(
            phred_scores,
            threshold=20
        )

        (
            trimmed_sequence,
            trimmed_quality_string,
            trimmed_scores
        ) = trim_low_quality_ends(
            sequence,
            quality_string,
            phred_scores,
            threshold=20
        )

        report = create_report(
            filename=filename,
            original_sequence_length=original_sequence_length,
            original_quality_length=original_quality_length,
            sequence=sequence,
            quality_string=quality_string,
            selected_encoding=selected_encoding,
            min_ascii=min_ascii,
            max_ascii=max_ascii,
            statistics=statistics,
            low_quality_positions=low_quality_positions,
            low_quality_regions=low_quality_regions,
            trimmed_sequence=trimmed_sequence,
            trimmed_quality_string=trimmed_quality_string,
            trimmed_scores=trimmed_scores,
            phred_scores=phred_scores
        )

        print("\n" + report)

        if selected_encoding == "Phred+64":
            convert_choice = input(
                "\nWould you like to convert the quality string "
                "to Phred+33? (y/n): "
            ).strip().lower()

            if convert_choice == "y":
                converted_quality = phred_to_quality_string(
                    phred_scores,
                    "Phred+33"
                )

                print(
                    "\nQuality string converted to Phred+33:"
                )
                print(converted_quality)

        save_choice = input(
            "\nWould you like to save the report? (y/n): "
        ).strip().lower()

        if save_choice == "y":
            output_filename = input(
                "Enter the output filename: "
            ).strip()

            save_report(
                report,
                output_filename
            )

    except (
        ValueError,
        FileNotFoundError,
        OSError
    ) as error:
        print(f"\nError: {error}")


if __name__ == "__main__":
    main()