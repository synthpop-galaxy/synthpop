from pydantic import BaseModel as _BaseModel
from pydantic import validator, Extra, root_validator
from typing import Literal, Any, Union, Callable
from pydantic.typing import AnyCallable


def model_validator(
        *,
        mode: Literal['wrap', 'before', 'after'],
        ) -> Union['AnyClassMethod', Callable[[AnyCallable], 'AnyClassMethod']]:
    return root_validator(pre=True)


class BaseModel(_BaseModel):
    @property
    def model_extra(self):
        attributes = {k: v for k, v in self.__dict__.items() if k not in self.__fields__}
        return attributes

    @classmethod
    def model_validate(cls, obj):
        return cls.parse_obj(obj)

    def model_dump_json(self, *args, **kwargs):
        return self.json(*args, **kwargs)
