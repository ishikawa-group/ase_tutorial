```
ase.io.read(filename: Union[str, PurePath, IO], index: Any = None, format: Optional[str] = None, parallel: bool = True, do_not_split_by_at_sign: bool = False, **kwargs)→ Union[Atoms, List[Atoms]][source]¶
Read Atoms object(s) from file.

filename: str or file
Name of the file to read from or a file descriptor.

index: int, slice or str
The last configuration will be returned by default. Examples:

index=0: first configuration

index=-2: second to last

index=':' or index=slice(None): all

index='-3:' or index=slice(-3, None): three last

index='::2' or index=slice(0, None, 2): even

index='1::2' or index=slice(1, None, 2): odd
```