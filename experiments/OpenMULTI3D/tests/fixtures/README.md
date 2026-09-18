# Test fixtures

Small example inputs that are read by the tests and reviewed like source.

Admission rules:

- A file belongs here when its content can be checked by reading it, or when it was
  produced by a recorded command whose provenance is written in the test that uses it.
- Keep it small. Hand-written text and a few kilobytes of binary are the intended size.
  A model atmosphere, a full line list or a run output directory belongs outside the
  repository; a test that needs one carries the `requires_input_data` marker.
- Record units, axis order and the wavelength medium wherever the format does not.
- Generated outputs stay out. Nothing here is written by a test run.
